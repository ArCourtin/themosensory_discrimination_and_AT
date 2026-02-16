# Script used to assess the model recovery results
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)

rm(list=ls())

#### Extract and aggregate data ####
model_data <-
  read_csv("recovery_analysis/absolute_model_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature
  )

model_data <-
  read_csv("recovery_analysis/relative_model_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+50
  ) %>%
  full_join(model_data)

model_data <-
  read_csv("recovery_analysis/mixed_model_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+100
  ) %>%
  full_join(model_data) %>% 
  filter(trial==1,participant==1,adapting_temperature_idx==1) %>% 
  pivot_longer(
    cols=c('mu_log_alpha','mu_log_beta','mu_hlogit_lambda','tau_log_alpha','tau_log_beta','tau_hlogit_lambda','mu_log_sigma','tau_log_sigma'),
    names_to='variable',
    values_to = 'truth'
    )


#### Loop through datasets - absolute coding ####
result_summary<-NULL
for(dataset in 1:50){
  result_summary <- readRDS(paste0("recovery_analysis/results/fits/summary_1_1_",dataset,".rds")) %>% 
    mutate(
      dataset=dataset,
      model='a',
      variable=c('mu_log_alpha','mu_log_beta','mu_hlogit_lambda','tau_log_alpha','tau_log_beta','tau_hlogit_lambda')
      ) %>% 
    bind_rows(result_summary)
}
for(dataset in 51:100){
  result_summary <- readRDS(paste0("recovery_analysis/results/fits/summary_2_2_",dataset,".rds")) %>% 
    mutate(
      dataset=dataset,
      model='r',
      variable=c('mu_log_alpha','mu_log_beta','mu_hlogit_lambda','tau_log_alpha','tau_log_beta','tau_hlogit_lambda')
    ) %>%     
    bind_rows(result_summary)
}
for(dataset in 101:150){
  result_summary <- readRDS(paste0("recovery_analysis/results/fits/summary_3_3_",dataset,".rds")) %>% 
    mutate(
      dataset=dataset,
      model='m',
      variable=c('mu_log_alpha','mu_log_sigma','mu_hlogit_lambda','tau_log_alpha','tau_log_sigma','tau_hlogit_lambda')
    ) %>%     
    bind_rows(result_summary)
}

pooled<-result_summary %>% 
  full_join(model_data) %>% 
  filter(!is.na(truth)) %>% 
  mutate(model=factor(model,c('a','r','m'),c('absolute','relative','mixed')))

pooled %>% 
  ggplot()+
  geom_pointrange(
    aes(
      x=truth,
      y=mean,
      ymin=q5,
      ymax=q95
        ),
    size=.1
    )+
  geom_abline()+
  facet_wrap(model~variable,scales = 'free',nrow = 3)+
  theme_classic()+
  labs(x='True value', y='Estimated value')



