# Script to used to iteratively fit the different models to the different simulated datasets
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(loo)
library(furrr)
library(future)
library(future.apply)


rm(list=ls())

Sys.setenv(
  PATH = "/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin",  
  CXX = "g++-9",
  CXX14 = "g++-9"
)

if (!dir.exists("recovery_analysis/sampling/c++_models")) {
  dir.create("recovery_analysis/sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/sampling/output")) {
  dir.create("recovery_analysis/sampling/output", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/sampling/draws")) {
  dir.create("recovery_analysis/sampling/draws", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/results/fits")) {
  dir.create("recovery_analysis/results/fits", recursive = TRUE)
}
if (!dir.exists("recovery_analysis/results/loo")) {
  dir.create("recovery_analysis/results/loo", recursive = TRUE)
}

#### Functions ####
fit_model <- function(iter_info) {
  file_name_out <- paste0("recovery_analysis/sampling/output/Output_rating_",iter_info$generative,'_',iter_info$fitted,'_',iter_info$dataset,".txt")
  
  con_out <- file(file_name_out, open = "wt")
  
  sink(con_out, append = TRUE)
  
  on.exit({
    sink(NULL)
    close(con_out)
  }, add = TRUE)
  
  sample_discrimination_data <-
    discrimination_data %>%
    filter(dataset == iter_info$dataset)
  sample_rating_data <-
    rating_data %>%
    filter(dataset == iter_info$dataset)
  
  data_list <- list(
    N  = nrow(sample_rating_data),
    P  = length(unique(sample_rating_data$participant)),
    is_cold = 0,
    recorded_baseline_temperature =
      sample_rating_data$recorded_baseline_temperature,
    
    absolute_target_temperature =
      sample_rating_data$absolute_target_temperature,
    
    absolute_adapting_temperature =
      sample_rating_data$absolute_adapting_temperature,
    
    rating =
      sample_rating_data$rating,
    
    participant =
      sample_rating_data$participant,
    
    # Discrimination block
    D = nrow(sample_discrimination_data),
    
    recorded_baseline_temperature_disc =
      sample_discrimination_data$recorded_baseline_temperature,
    
    absolute_target_temperature_disc =
      sample_discrimination_data$absolute_target_temperature,
    
    absolute_adapting_temperature_disc =
      sample_discrimination_data$absolute_adapting_temperature,
    
    choice_accuracy =
      sample_discrimination_data$choice_accuracy,
    
    participant_disc =
      sample_discrimination_data$participant
  )
  
  
  mod <- compiled_models[[iter_info$model]]
  
  seed   <- 12345
  mt     <- 12
  ad     <- 0.99
  n_iter <- 2000
  chains <- 4
  
  fit <- mod$sample(
    data = data_list,
    seed = seed,
    chains = chains,
    parallel_chains = chains,
    iter_warmup = n_iter,
    iter_sampling = n_iter,
    save_warmup = FALSE,
    max_treedepth = mt,
    adapt_delta = ad,
    refresh = n_iter / 10
  )
  gc()
  
  d<-fit$diagnostic_summary()
  print(d)
  saveRDS(d,paste0('recovery_analysis/results/fits/diagnostics_rating_',iter_info$generative,'_',iter_info$fitted,'_',iter_info$dataset,'.rds'))
  gc()
  
  s<-fit$summary(c('mu','tau'))
  print(s)
  saveRDS(s,paste0('recovery_analysis/results/fits/summary_rating_',iter_info$generative,'_',iter_info$fitted,'_',iter_info$dataset,'.rds'))
  gc()
  
  loo<-fit$loo(cores=chains)
  saveRDS(loo,paste0('recovery_analysis/results/loo/loo_rating_',iter_info$generative,'_',iter_info$fitted,'_',iter_info$dataset,'.rds'))
  gc()
}
#### Compile models ####
models<-
  c(
    'rating_absolute_coding.stan',
    'rating_relative_coding.stan',
    'rating_mixed_coding.stan',
    'rating_non_mechanistic.stan'
  )

plan(multisession, workers = 4)

compiled_models <- future_lapply(models, function(m) {
  cmdstan_model(
    stan_file = file.path("stan_models", m),
    dir = "recovery_analysis/sampling/c++_models",
    stanc_options = list("O1"),
    compile_model_methods = F,
    force_recompile = F
  )
})

names(compiled_models) <- models

#### Extract and aggregate data ####
discrimination_data <-
  read_csv("recovery_analysis/simulated_data/absolute_model_discrimination_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature
  )

discrimination_data <-
  read_csv("recovery_analysis/simulated_data/relative_model_discrimination_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+50
  ) %>%
  full_join(discrimination_data)

discrimination_data <-
  read_csv("recovery_analysis/simulated_data/mixed_model_discrimination_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+100
  ) %>%
  full_join(discrimination_data)
discrimination_data<-discrimination_data

rating_data <-
  read_csv("recovery_analysis/simulated_data/absolute_model_rating_data.csv") %>%
  mutate(
    relative_adapting_temperature =
    absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature
    )

rating_data <-
  read_csv("recovery_analysis/simulated_data/relative_model_rating_data.csv") %>%
  mutate(
    relative_adapting_temperature =
    absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+50
    ) %>%
  full_join(rating_data)

rating_data <-
  read_csv("recovery_analysis/simulated_data/mixed_model_rating_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+100
    ) %>%
  full_join(rating_data)

#### Prepare lists for fitting runs ##############
mod_comb<-expand_grid(dataset=1:150,fitted=1:4)
iter_info<-list()
for(m in 1:dim(mod_comb)[1]){
  iter_info[[m]]<-
    list(
      generative=ceiling(mod_comb$dataset[m]/50),
      fitted=mod_comb$fitted[m],
      model=models[mod_comb$fitted[m]],
      dataset=mod_comb$dataset[m]
    )
}

plan(multisession, workers = 4)

# results<-future_map(
#   iter_info,
#   ~fit_model(.x),
#   .options = furrr_options(
#     seed = T,
#     scheduling = F,
#     stdout = F,
#     conditions = character()),
#   .progress = F)
