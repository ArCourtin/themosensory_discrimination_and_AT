# Script to used to iteratively fit the different models to the different simulated datasets
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(cmdstanr)
library(tidyverse)
library(loo)

rm(list=ls())

directory=getwd()

Sys.setenv(
  PATH = "/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin",  
  CXX = "g++-9",
  CXX14 = "g++-9"
)

if (!dir.exists("sampling/c++_models")) {
  dir.create("sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("sampling/output")) {
  dir.create("sampling/output", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/fits")) {
  dir.create("recovery_analysis/fits", recursive = TRUE)
}

#### Import and reformat data ####
absolute_model_data <- 
  read_csv("recovery_analysis/absolute_model_data.csv") %>% 
  mutate(
    relative_adapting_temperature=absolute_adapting_temperature-recorded_baseline_temperature,
    adapting_temperature_idx=3+relative_adapting_temperature
    )
relative_model_data <- 
  read_csv("recovery_analysis/relative_model_data.csv") %>% 
  mutate(
    relative_adapting_temperature=absolute_adapting_temperature-recorded_baseline_temperature,
    adapting_temperature_idx=3+relative_adapting_temperature
  )
mixed_model_data <- 
  read_csv("recovery_analysis/mixed_model_data.csv") %>% 
  mutate(
    relative_adapting_temperature=absolute_adapting_temperature-recorded_baseline_temperature,
    adapting_temperature_idx=3+relative_adapting_temperature
  )

#subset data
sample_data<-
  absolute_model_data %>% 
  filter(dataset==d)

#format data
data_list<-
  list(
    N=dim(sample_data)[1],
    P=length(unique(sample_data$participant)),
    recorded_baseline_temperature=sample_data$recorded_baseline_temperature,
    absolute_target_temperature=sample_data$absolute_target_temperature,
    absolute_adapting_temperature=sample_data$absolute_adapting_temperature,
    choice_accuracy=sample_data$choice_accuracy,
    participant=sample_data$participant,
    adapting_temperature_idx=sample_data$adapting_temperature_idx
  )

#### prepare models ####
mod_A<-cmdstan_model(stan_file='stan_models/discrimination_absolute_coding.stan',
                   dir="sampling/c++_models",
                   stanc_options = list("O1"),
                   force_recompile=T,
                   compile_model_methods=T
)
mod_R<-cmdstan_model(stan_file='stan_models/discrimination_relative_coding.stan',
                    dir="sampling/c++_models",
                    stanc_options = list("O1"),
                    force_recompile=T,
                    compile_model_methods=T
)
mod_M<-cmdstan_model(stan_file='stan_models/discrimination_mixed_coding.stan',
                   dir="sampling/c++_models",
                   stanc_options = list("O1"),
                   force_recompile=T,
                   compile_model_methods=T
)
mod_N<-cmdstan_model(stan_file='stan_models/discrimination_non_mechanistic.stan',
                   dir="sampling/c++_models",
                   stanc_options = list("O1"),
                   force_recompile=T,
                   compile_model_methods=T
)

seed=12345
mt=12
ad=.99
n_iter=1000
chains=4

#absolute coding model
pathfinder_fit<-
  mod_A$pathfinder(
    data=data_list,
    psis_resample = F,
    calculate_lp = F,
    seed = seed,
    refresh=500
  )

fit_A<-mod_A$sample(
  data = data_list,
  seed=seed,
  chains=chains,
  parallel_chains = chains,
  iter_warmup =n_iter,
  iter_sampling= n_iter,
  save_warmup=F,
  max_treedepth = mt,
  adapt_delta = ad,
  refresh = n_iter/10,
  init=pathfinder_fit
)
d_A<-fit_A$diagnostic_summary() 
s_A<-fit_A$summary(c('mu','tau')) 
loo_A<-fit_A$loo(moment_match = T)

#relative coding model
pathfinder_fit<-
  mod_R$pathfinder(
    data=data_list,
    psis_resample = F,
    calculate_lp = F,
    seed = seed,
    refresh=500
  )

fit_R<-mod_R$sample(
  data = data_list,
  seed=seed,
  chains=chains,
  parallel_chains = chains,
  iter_warmup =n_iter,
  iter_sampling= n_iter,
  save_warmup=F,
  max_treedepth = mt,
  adapt_delta = ad,
  refresh = n_iter/10,
  init=pathfinder_fit
)
d_R<-fit_R$diagnostic_summary() 
s_R<-fit_R$summary(c('mu','tau')) 
loo_R<-fit_R$loo(moment_match = T)

#mixed coding model
pathfinder_fit<-
  mod_M$pathfinder(
    data=data_list,
    psis_resample = F,
    calculate_lp = F,
    seed = seed,
    refresh=500
  )

fit_M<-mod_M$sample(
  data = data_list,
  seed=seed,
  chains=chains,
  parallel_chains = chains,
  iter_warmup =n_iter,
  iter_sampling= n_iter,
  save_warmup=F,
  max_treedepth = mt,
  adapt_delta = ad,
  refresh = n_iter/10,
  init=pathfinder_fit
)
d_M<-fit_M$diagnostic_summary() 
s_M<-fit_M$summary(c('mu','tau')) 
loo_M<-fit_M$loo(moment_match = T)

#non-mechansitic model
pathfinder_fit<-
  mod_N$pathfinder(
    data=data_list,
    psis_resample = F,
    calculate_lp = F,
    seed = seed,
    refresh=500
  )

fit_N<-mod_N$sample(
  data = data_list,
  seed=seed,
  chains=chains,
  parallel_chains = chains,
  iter_warmup =n_iter,
  iter_sampling= n_iter,
  save_warmup=F,
  max_treedepth = mt,
  adapt_delta = ad,
  refresh = n_iter/10,
  init=pathfinder_fit
)
d_N<-fit_N$diagnostic_summary() 
s_N<-fit_N$summary(c('mu','tau')) 
loo_N<-fit_N$loo(moment_match = T)

#compare models
mod_comp<-loo_compare(
  list(
    absolute=loo_A,
    relative=loo_R,
    mixed=loo_M,
    non_mech=loo_N
  )
) %>% 
  as_tibble() %>% 
  mutate(
    z=elpd_diff/se_diff,
    p=pnorm(z)
    )


plan(multisession, workers = 4)
results<-future_map(
  parameters_fit,
  ~possfit_model(.x),
  .options = furrr_options(
    seed = T,
    scheduling = F,
    stdout = F,
    conditions = character()),
  .progress = F)










