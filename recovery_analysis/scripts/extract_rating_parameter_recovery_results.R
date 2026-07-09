# Script used to assess the parameter recovery results for the rating models
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)

rm(list=ls())

#### Extract and aggregate data ####
model_data <-
  read_csv("recovery_analysis/simulated_data/absolute_model_rating_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature
  )

model_data <-
  read_csv("recovery_analysis/simulated_data/relative_model_rating_data.csv") %>%
  mutate(
    relative_adapting_temperature =
      absolute_adapting_temperature - recorded_baseline_temperature,
    adapting_temperature_idx = 3 + relative_adapting_temperature,
    dataset=dataset+50
  ) %>%
  full_join(model_data) %>% 
  filter(trial==1,participant==1,adapting_temperature_idx==1) %>% 
  pivot_longer(
    cols=c('mu_intercept','mu_log_slope','mu_log_lb','mu_log_ub','mu_log_eta','tau_intercept','tau_log_slope','tau_log_lb','tau_log_ub','tau_log_eta'),
    names_to='variable',
    values_to = 'truth'
    )


#### Loop through datasets - absolute coding ####
result_summary<-NULL
for(dataset in 1:50){
  result_summary <- readRDS(paste0("recovery_analysis/results/fits/summary_rating_1_1_",dataset,".rds")) %>% 
    mutate(
      dataset=dataset,
      model='a',
      variable=c('mu_intercept','mu_log_slope','mu_log_lb','mu_log_ub','mu_log_eta','tau_intercept','tau_log_slope','tau_log_lb','tau_log_ub','tau_log_eta')
      ) %>% 
    bind_rows(result_summary)
}
for(dataset in 51:100){
  result_summary <- readRDS(paste0("recovery_analysis/results/fits/summary_rating_2_2_",dataset,".rds")) %>% 
    mutate(
      dataset=dataset,
      model='r',
      variable=c('mu_intercept','mu_log_slope','mu_log_lb','mu_log_ub','mu_log_eta','tau_intercept','tau_log_slope','tau_log_lb','tau_log_ub','tau_log_eta')
    ) %>%     
    bind_rows(result_summary)
}

pooled<-
  result_summary %>% 
  full_join(model_data) %>% 
  filter(!is.na(truth)) %>% 
  mutate(
    model=factor(model,c('a','r'),c('absolute','relative')),
    variable=factor(variable,c('mu_intercept','mu_log_slope','mu_log_lb','mu_log_ub','mu_log_eta','tau_intercept','tau_log_slope','tau_log_lb','tau_log_ub','tau_log_eta'))
    )

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
  facet_wrap(model~variable,scales = 'free',ncol = 5)+
  theme_classic()+
  labs(x='True value', y='Estimated value')


ggsave('recovery_analysis/figures/PR_R.png',units = 'cm',width = 30,height = 20)

