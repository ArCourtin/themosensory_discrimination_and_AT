# Script to used to iteratively fit the different models to the different simulated datasets
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(rslurm)
library(loo)

rm(list=ls())

wd<-getwd()

if (!dir.exists("sampling/c++_models")) {
  dir.create("sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("sampling/output")) {
  dir.create("sampling/output", recursive = TRUE)
}
if (!dir.exists("results/fits")) {
  dir.create("results/fits", recursive = TRUE)
}
if (!dir.exists("results/loo")) {
  dir.create("results/loo", recursive = TRUE)
}

#### Functions ####
fit_model <- function(iter_info) {

  base_dir=iter_info$wd
  print(iter_info$model_path)
  mod <-   cmdstan_model(
    stan_file = iter_info$model_path,
    dir = file.path(base_dir,"sampling","c++_models"),
    stanc_options = list("O1"),
    force_recompile=T,
    compile_model_methods=T
    )
  seed=1234
  pathfinder_fit<-
    mod$pathfinder(
      data=iter_info$data,
      psis_resample = F,
      calculate_lp = F,
      seed = seed,
      refresh=500
    )
  
  fit <- mod$sample(
    data = iter_info$data_list,
    seed = seed,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    max_treedepth = 12,
    adapt_delta = 0.95,
    init=pathfinder_fit,
    refresh = 200
  )
  fit$diagnostic_summary() %>% print()
  fit$summary(c('mu','tau')) %>% print(n=50)
  fit$save_object(file.path(base_dir, "results","fits",paste0(iter_info$fitted,'_',iter_info$task,".rds")))
  
  loo <- fit$loo(cores = 4,moment_match = T)
  print(loo)
  saveRDS(loo, file.path(base_dir, "results","loo",paste0(iter_info$fitted,'_',iter_info$task,".rds")))
}

#### Compile models ####

model_paths<-
  c(
    file.path(wd,"stan_models","discrimination_absolute_coding.stan"),
    file.path(wd,"stan_models","discrimination_relative_coding.stan"),
    file.path(wd,"stan_models","discrimination_non_mechanistic.stan")
  )

compiled_models <- lapply(model_paths, function(p) {
  cmdstan_model(
    stan_file = p,
    dir = file.path(wd, "sampling","c++_models"),
    stanc_options = list("O1")
  )
})

names(compiled_models) <- model_paths

#### Extract data ####
data <-
  read_csv("data/d_at_2ifc.csv")%>% 
  mutate(
    relative_adapting_temperature =
    round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature
    )
participant<-unique(data$participant)
P<-length(participant)
for(pdx in 1:P){
  data$participant[data$participant==participant[pdx]]<-pdx
}
#### Prepare lists for fitting runs ##############
iter_info=list()
for(t in 0:1){
  sample_data<-data %>% filter(task==t)
  for(m in 1:3){
    data_list <- list(
      N = nrow(sample_data),
      P = length(unique(sample_data$participant)),
      recorded_baseline_temperature = sample_data$baseline,
      absolute_target_temperature   = sample_data$adapting+t*sample_data$temperature+(t-1)*sample_data$temperature,
      absolute_adapting_temperature = sample_data$adapting,
      choice_accuracy               = sample_data$accuracy,
      participant                   = sample_data$participant,
      adapting_temperature_idx      = sample_data$adapting_temperature_idx,
      is_cold = t==0
    )
    
    iter_info[[m+t*3]]<-
      list(
        wd=wd,
        fitted=m,
        task=t+1,
        model_path = model_paths[m],
        data_list=data_list
      )
  }
}


#### Launch slurm jobs ####
dir.create("slurm", showWarnings = FALSE)
setwd("slurm")

jobs <- vector("list", length(iter_info))

for (k in seq_along(iter_info)) {
  jobs[[k]] <- slurm_map(
    x = list(iter_info[[k]]),
    f = fit_model,
    jobname = paste0("stan_", k),
    cpus_per_node = 4,
    nodes = 1,
    slurm_options = list(time = "24:00:00", mem = "8G"),
    pkgs = c("tidyverse", "cmdstanr", "loo")
  )
}

setwd(wd)
