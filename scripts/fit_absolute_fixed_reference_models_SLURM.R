# Companion fit to the main absolute-coding models (discrimination_absolute_coding.stan,
# rating_absolute_coding.stan), testing whether the absolute reference is better described as
# each participant's own measured pre-session baseline (model 1) or a fixed value common to
# everyone (model 2, fit here). Reuses both stan files completely unchanged: the only difference
# is that recorded_baseline_temperature is set to a constant 32 for every row instead of each
# participant's measured baseline, so the model itself is identical, just fed a different
# reference. Saved under model index 2 (results/{fits,loo}/{discrimination,rating}_2_{task}.rds),
# alongside models 1 (absolute, personal baseline), 3 (relative), and 4 (non-mechanistic) from
# fit_discrimination_models_SLURM.R / fit_rating_models_SLURM.R, so all four can be compared
# directly via model_comparison.R.
#
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Written with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(rslurm)
library(loo)

rm(list=ls())

wd<-getwd()
fixed_reference<-32

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
  fit$save_object(file.path(base_dir, "results","fits",paste0(iter_info$domain,"_2_",iter_info$task,".rds")))

  loo <- fit$loo(cores = 4,moment_match = T)
  print(loo)
  saveRDS(loo, file.path(base_dir, "results","loo",paste0(iter_info$domain,"_2_",iter_info$task,".rds")))
}

#### Compile models ####

model_paths<-c(
  discrimination=file.path(wd,"stan_models","discrimination_absolute_coding.stan"),
  rating=file.path(wd,"stan_models","rating_absolute_coding.stan")
)

compiled_models <- lapply(model_paths, function(p) {
  cmdstan_model(
    stan_file = p,
    dir = file.path(wd, "sampling","c++_models"),
    stanc_options = list("O1")
  )
})

#### Extract data (discrimination) ####
# Response coding: `active_interval` gives the 2IFC interval that held the deviating stimulus
# (1 = first, 2 = second). The second-interval choice is reconstructed from accuracy, i.e.
# correct & active in interval 2 -> chose 2, correct & active in interval 1 -> chose 1, etc.
# baseline_flag==1 marks trials where the probe's pre-trial baseline temperature drifted more
# than 0.2C from its mean, i.e. the assumed 32C reference did not actually hold for that trial;
# these are excluded before fitting.
discrimination_data <-
  read_csv("data/d_at_2ifc_af.csv")%>%
  filter(baseline_flag == 0, deviation_flag==0) %>%
  mutate(
    relative_adapting_temperature =
    round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    interval_sign = if_else(active_interval == 2, 1, -1),
    chose_second  = as.integer(if_else(active_interval == 2, accuracy, 1L - accuracy))
    )
participant<-unique(discrimination_data$participant)
for(pdx in seq_along(participant)){
  discrimination_data$participant[discrimination_data$participant==participant[pdx]]<-pdx
}

#### Extract data (rating) ####
rating_data <-
  read_csv("data/d_at_ratings_af.csv") %>%
  filter(baseline_flag == 0, deviation_flag==0, confirmed == 1,!is.na(recorded_temperature)) %>%
  mutate(
    relative_adapting_temperature = round(adapting - baseline),
    adapting_temperature_idx = 3 + relative_adapting_temperature
  )
participant<-unique(rating_data$participant)
for(pdx in seq_along(participant)){
  rating_data$participant[rating_data$participant==participant[pdx]]<-pdx
}

#### Prepare lists for fitting runs (absolute model only, both domains, both tasks) ####
iter_info=list()
idx<-1
for(t in 0:1){
  d_sample<-discrimination_data %>% filter(task==t)
  iter_info[[idx]]<-list(
    wd=wd,
    domain="discrimination",
    task=t+1,
    model_path=model_paths[["discrimination"]],
    data_list=list(
      N = nrow(d_sample),
      P = length(unique(d_sample$participant)),
      recorded_baseline_temperature = rep(fixed_reference, nrow(d_sample)),
      absolute_target_temperature   = round(d_sample$recorded_temperature,1),
      absolute_adapting_temperature = d_sample$adapting,
      interval_sign                 = d_sample$interval_sign,
      chose_second                  = d_sample$chose_second,
      participant                   = d_sample$participant,
      adapting_temperature_idx      = d_sample$adapting_temperature_idx,
      is_cold = t==0
    )
  )
  idx<-idx+1

  r_sample<-rating_data %>% filter(task==t)
  iter_info[[idx]]<-list(
    wd=wd,
    domain="rating",
    task=t+1,
    model_path=model_paths[["rating"]],
    data_list=list(
      N = nrow(r_sample),
      P = length(unique(r_sample$participant)),
      recorded_baseline_temperature = rep(fixed_reference, nrow(r_sample)),
      absolute_target_temperature   = round(r_sample$recorded_temperature,1),
      rating                        = r_sample$rating/100,
      participant                   = r_sample$participant,
      is_cold = t==0
    )
  )
  idx<-idx+1
}

#### Launch slurm jobs ####
dir.create("slurm", showWarnings = FALSE)
setwd("slurm")

jobs <- vector("list", length(iter_info))

for (k in seq_along(iter_info)) {
  jobs[[k]] <- slurm_map(
    x = list(iter_info[[k]]),
    f = fit_model,
    jobname = paste0("stan_fixedref_", k),
    cpus_per_node = 4,
    nodes = 1,
    slurm_options = list(time = "24:00:00", mem = "8G"),
    pkgs = c("tidyverse", "cmdstanr", "loo")
  )
}

setwd(wd)
