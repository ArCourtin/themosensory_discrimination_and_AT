# Script to used to iteratively fit the different models to the different simulated datasets
# Model indices used for saved filenames: 1=absolute (personal-baseline reference), 3=relative,
# 4=non-mechanistic. Index 2 (absolute, fixed common reference) is fit separately by
# fit_absolute_fixed_reference_models_SLURM.R, reusing the same absolute stan model unchanged.
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic).

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
      data=iter_info$data_list,
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
    iter_warmup = 2000,
    iter_sampling = 2000,
    max_treedepth = 12,
    adapt_delta = 0.99,
    init=pathfinder_fit,
    refresh = 200
  )
  fit$diagnostic_summary() %>% print()
  fit$summary(c('mu','tau')) %>% print(n=50)
  fit$save_object(file.path(base_dir, "results","fits",paste0("discrimination_",iter_info$fitted,'_',iter_info$task,".rds")))
  
  loo <- fit$loo(cores = 4,moment_match = T)
  print(loo)
  saveRDS(loo, file.path(base_dir, "results","loo",paste0("discrimination_",iter_info$fitted,'_',iter_info$task,".rds")))
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

#### Extract data ####
# Response coding: `active_interval` gives the 2IFC interval that held the deviating stimulus
# (1 = first, 2 = second). The second-interval choice is reconstructed from accuracy, i.e.
# correct & active in interval 2 -> chose 2, correct & active in interval 1 -> chose 1, etc.
# baseline_flag==1 marks trials where the probe's pre-trial baseline temperature drifted more
# than 0.2C from its mean, i.e. the assumed 32C reference did not actually hold for that trial;
# these are excluded before fitting.
data <-
  read_csv("data/d_at_2ifc.csv")%>%
  filter(baseline_flag == 0) %>%
  mutate(
    relative_adapting_temperature =
    round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    interval_sign = if_else(active_interval == 2, 1, -1),
    chose_second  = as.integer(if_else(active_interval == 2, accuracy, 1L - accuracy))
    )
participant<-unique(data$participant)
P<-length(participant)
for(pdx in 1:P){
  data$participant[data$participant==participant[pdx]]<-pdx
}
#### Prepare lists for fitting runs ##############
# model_paths order is (absolute, relative, non-mechanistic); model_save_idx maps that compile
# order onto the saved-filename index, leaving slot 2 free for the fixed-reference absolute
# variant fit by fit_absolute_fixed_reference_models_SLURM.R.
model_save_idx<-c(1,3,4)
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
      interval_sign                 = sample_data$interval_sign,
      chose_second                  = sample_data$chose_second,
      participant                   = sample_data$participant,
      adapting_temperature_idx      = sample_data$adapting_temperature_idx,
      is_cold = t==0
    )

    iter_info[[m+t*3]]<-
      list(
        wd=wd,
        fitted=model_save_idx[m],
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
    jobname = paste0("stan_d_", k),
    cpus_per_node = 4,
    nodes = 1,
    slurm_options = list(time = "24:00:00", mem = "8G"),
    pkgs = c("tidyverse", "cmdstanr", "loo")
  )
}

setwd(wd)
