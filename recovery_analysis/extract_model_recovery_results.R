# Script used to assess the model recovery results
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(loo)

rm(list=ls())

#### Loop through iterations to extract informations - full fitted model space ####
models<-c('absolute','relative','mixed','non-mechanistic')
mod_comp_res<-NULL
for(generative in 1:3){
  for(dataset in 1:50){
    div_a<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_1_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    div_r<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_2_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    div_m<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_3_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    div_n<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_4_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    
    loo_a<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_1_",dataset+(generative-1)*50,".rds"))
    loo_r<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_2_",dataset+(generative-1)*50,".rds"))
    loo_m<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_3_",dataset+(generative-1)*50,".rds"))
    loo_n<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_4_",dataset+(generative-1)*50,".rds"))
    
    k_a<-sum(loo_a$diagnostics[[1]]>.7)
    k_r<-sum(loo_r$diagnostics[[1]]>.7)
    k_m<-sum(loo_m$diagnostics[[1]]>.7)
    k_n<-sum(loo_n$diagnostics[[1]]>.7)
    
    comp<-loo_compare(
      list(
        absolute=loo_a,
        relative=loo_r,
        mixed=loo_m,
        `non-mechanistic`=loo_n
        )
      )
    rn<-row.names(comp)
    comp<-
      comp%>% 
      as_tibble() %>% 
      mutate(
        fitted=rn,
        z=elpd_diff/se_diff,
        p=pnorm(z),
        won=is.na(z),
        ) %>% 
      arrange(fitted)
    comp_sig<-
      comp %>% 
      mutate(sig=p<0.05) %>% 
      summarise(sig=3==sum(sig,na.rm = T))
    
    comp$won_sig<-comp_sig$sig
    comp$div<-c(div_a,div_m,div_n,div_r)
    comp$k<-c(k_a,k_m,k_n,k_r)
    comp$generative<-models[generative]
    comp$dataset<-dataset+(generative-1)*50
    
    mod_comp_res<-bind_rows(mod_comp_res,comp)
  }
}

results<-
  mod_comp_res %>% 
  mutate(
    generative=factor(generative,c('absolute','relative','mixed')),
    fitted=factor(fitted,c('absolute','relative','mixed','non-mechanistic'))
    ) %>% 
  group_by(generative,fitted) %>% 
  summarise(n=sum(won)) %>% 
  mutate(
    ratio=n/50,
    lb=qbeta(.025,1+n,1+50-n),
    ub=qbeta(.975,1+n,1+50-n)
  )

results %>% 
  ggplot(aes(x=generative,y=fitted,fill=ratio))+
  geom_tile()+
  geom_text(
    aes(
      label = sprintf(
        "%.2f\n[%.2f, %.2f]",
        ratio, lb, ub
      )
    ),
    size = 3
  ) +
  scale_fill_viridis_c(limits=c(0,1)) +
  theme_minimal()


#### Loop through iterations to extract informations - restricted fitted model space ####
models<-c('absolute','relative','mixed','non-mechanistic')
mod_comp_res<-NULL
for(generative in 1:3){
  for(dataset in 1:50){
    div_a<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_1_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    div_r<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_2_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    div_m<-
      readRDS(paste0("recovery_analysis/results/fits/diagnostics_",generative,"_3_",dataset+(generative-1)*50,".rds"))[[1]] %>% 
      sum()
    
    loo_a<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_1_",dataset+(generative-1)*50,".rds"))
    loo_r<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_2_",dataset+(generative-1)*50,".rds"))
    loo_m<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_3_",dataset+(generative-1)*50,".rds"))

    k_a<-sum(loo_a$diagnostics[[1]]>.7)
    k_r<-sum(loo_r$diagnostics[[1]]>.7)
    k_m<-sum(loo_m$diagnostics[[1]]>.7)

    comp<-loo_compare(
      list(
        absolute=loo_a,
        relative=loo_r,
        mixed=loo_m
      )
    )
    rn<-row.names(comp)
    comp<-
      comp%>% 
      as_tibble() %>% 
      mutate(
        fitted=rn,
        z=elpd_diff/se_diff,
        p=pnorm(z),
        won=is.na(z),
      ) %>% 
      arrange(fitted)
    comp_sig<-
      comp %>% 
      mutate(sig=p<0.05) %>% 
      summarise(sig=2==sum(sig,na.rm = T))
    
    comp$won_sig<-comp_sig$sig
    comp$div<-c(div_a,div_m,div_r)
    comp$k<-c(k_a,k_m,k_r)
    comp$generative<-models[generative]
    comp$dataset<-dataset+(generative-1)*50
    
    mod_comp_res<-bind_rows(mod_comp_res,comp)
  }
}

results<-
  mod_comp_res %>% 
  mutate(
    generative=factor(generative,c('absolute','relative','mixed')),
    fitted=factor(fitted,c('absolute','relative','mixed'))
  ) %>% 
  group_by(generative,fitted) %>% 
  summarise(n=sum(won)) %>% 
  mutate(
    ratio=n/50,
    lb=qbeta(.025,1+n,1+50-n),
    ub=qbeta(.975,1+n,1+50-n)
  )

results %>% 
  ggplot(aes(x=generative,y=fitted,fill=ratio))+
  geom_tile()+
  geom_text(
    aes(
      label = sprintf(
        "%.2f\n[%.2f, %.2f]",
        ratio, lb, ub
      )
    ),
    size = 3
  ) +
  scale_fill_viridis_c(limits=c(0,1)) +
  theme_minimal()

#### Diagnose the confoudn for mixed ####
model_data <-
  read_csv("recovery_analysis/mixed_model_data.csv") %>%
  group_by(dataset) %>% 
  summarise(
    mu_log_alpha=mean(mu_log_alpha),
    mu_log_sigma=mean(mu_log_sigma),
    mu_hlogit_lambda=mean(mu_hlogit_lambda),
    tau_log_alpha=mean(tau_log_alpha),
    tau_log_sigma=mean(tau_log_sigma),
    tau_hlogit_lambda=mean(tau_hlogit_lambda)
  )
for(generative in 3){
  for(dataset in 1:50){
    loo_r<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_2_",dataset+(generative-1)*50,".rds"))
    loo_m<-readRDS(paste0("recovery_analysis/results/loo/loo_",generative,"_3_",dataset+(generative-1)*50,".rds"))

    model_data$elpd_diff[dataset]<-loo_m$estimates[1]-loo_r$estimates[1]
 }
}

model_data %>% 
  pivot_longer(cols=!c(dataset,elpd_diff)) %>% 
  ggplot()+
  geom_point(aes(x=value,y=elpd_diff))+
  facet_wrap(name~.,scales = 'free')+
  theme_classic()
