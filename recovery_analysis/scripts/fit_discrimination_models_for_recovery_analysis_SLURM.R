# Script to used to iteratively fit the different models to the different simulated datasets
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic). 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(rslurm)

rm(list=ls())

wd<-getwd()

# Number of simulated datasets per generative model (must match the MATLAB simulation).
# Datasets 1:n_datasets = absolute, (n_datasets+1):(2*n_datasets) = relative.
n_datasets <- 30

if (!dir.exists("recovery_analysis/sampling/c++_models")) {
  dir.create("recovery_analysis/sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/sampling/output")) {
  dir.create("recovery_analysis/sampling/output", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/results/fits")) {
  dir.create("recovery_analysis/results/fits", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/results/loo")) {
  dir.create("recovery_analysis/results/loo", recursive = TRUE)
}

#### Functions ####
fit_model <- function(iter_info) {

  base_dir=iter_info$wd
  
  mod <-   cmdstan_model(
    stan_file = iter_info$model_path,
    dir = file.path(base_dir, "recovery_analysis","sampling","c++_models"),
    stanc_options = list("O1")
  )
  
  seed=12345
  # Pathfinder can occasionally fail to initialize on a given simulated dataset; fall back to
  # cmdstanr's default random init (init=NULL) so a single failure does not lose that job.
  pathfinder_fit<-
    tryCatch(
      mod$pathfinder(
        data=iter_info$data_list,
        psis_resample = F,
        calculate_lp = F,
        seed = seed,
        refresh=500
      ),
      error = function(e) NULL
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
    refresh = 400
  )
  
  d <- fit$diagnostic_summary()
  saveRDS(d, file.path(base_dir, "recovery_analysis","results","fits",paste0("diagnostics_",iter_info$generative,"_",iter_info$fitted,"_",iter_info$dataset,".rds")))
  
  vars <- intersect(c("mu","tau"), fit$metadata()$stan_variables)
  s <- fit$summary(vars)
  saveRDS(s, file.path(base_dir, "recovery_analysis","results","fits",paste0("summary_",iter_info$generative,"_",iter_info$fitted,"_",iter_info$dataset,".rds")))
  
  loo <- fit$loo(cores = 4)
  saveRDS(loo, file.path(base_dir, "recovery_analysis","results","loo",paste0(iter_info$generative,"_",iter_info$fitted,"_",iter_info$dataset,".rds")))
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
    dir = file.path(wd, "recovery_analysis","sampling","c++_models"),
    stanc_options = list("O1")
  )
})

#### Extract and aggregate data ####
model_data <-
  read_csv("recovery_analysis/simulated_data/absolute_model_discrimination_data.csv") %>%
  mutate(
    relative_adapting_temperature =
    absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    interval_sign = if_else(active_interval == 2, 1, -1),
    chose_second  = as.integer(if_else(active_interval == 2, choice_accuracy, 1L - choice_accuracy))
    )

model_data <-
  read_csv("recovery_analysis/simulated_data/relative_model_discrimination_data.csv") %>%
  mutate(
    relative_adapting_temperature =
    absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    interval_sign = if_else(active_interval == 2, 1, -1),
    chose_second  = as.integer(if_else(active_interval == 2, choice_accuracy, 1L - choice_accuracy)),
    dataset=dataset+n_datasets
    ) %>%
  bind_rows(model_data)

#### Prepare lists for fitting runs ##############
mod_comb<-expand_grid(dataset=1:(2*n_datasets),fitted=1:3)
iter_info<-list()
for(m in 1:dim(mod_comb)[1]){
  sample_data <- model_data %>%
    filter(dataset == mod_comb$dataset[m])
  
  data_list <- list(
    N = nrow(sample_data),
    P = length(unique(sample_data$participant)),
    recorded_baseline_temperature = sample_data$recorded_baseline_temperature,
    absolute_target_temperature   = sample_data$absolute_target_temperature,
    absolute_adapting_temperature = sample_data$absolute_adapting_temperature,
    interval_sign                 = sample_data$interval_sign,
    chose_second                  = sample_data$chose_second,
    participant                   = sample_data$participant,
    adapting_temperature_idx      = sample_data$adapting_temperature_idx,
    is_cold = 0
  )
  
  iter_info[[m]]<-
    list(
      wd=wd,
      generative=ceiling(mod_comb$dataset[m]/n_datasets),
      fitted=mod_comb$fitted[m],
      dataset=mod_comb$dataset[m],
      model_path = model_paths[[mod_comb$fitted[m]]],
      data_list=data_list
    )
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
