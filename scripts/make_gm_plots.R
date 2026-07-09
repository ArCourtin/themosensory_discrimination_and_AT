# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)

rm(list=ls())
inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

#### Load fitted models ####
fit_abs <- readRDS("results/fits/1_1.rds")
fit_rel <- readRDS("results/fits/2_1.rds")
fit_nm <- readRDS("results/fits/3_1.rds")

#Extract posterior and reconstruct PF
pf_abs<-
  fit_abs$draws('mu',format = 'df') %>%
  mutate(
    rho=`mu[1]`,
    alpha=exp(`mu[2]`),
    beta=exp(`mu[3]`),
    lambda=.5*inv_logit(`mu[4]`)
  ) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,temperature=seq(28,35,.1),.draw=1:4000)) %>% 
  mutate(
    catt=baseline-temperature,
    caat=baseline-adapting,
    
    cs=catt-rho,
    sr=cs*inv_logit(100*cs),
    asr=sr*inv_logit(100*(sr-(caat+alpha-rho))),
    
    p=lambda+(1-2*lambda)*pnorm(beta*asr)
    ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='ab')

pf_abs %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

pf_rel<-
  fit_rel$draws('mu',format = 'df') %>%
  mutate(
    alpha=exp(`mu[1]`),
    beta=exp(`mu[2]`),
    lambda=.5*inv_logit(`mu[3]`)
  ) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,temperature=seq(28,35,.1),.draw=1:4000)) %>% 
  mutate(
    dfat=adapting-temperature,

    cs=dfat-alpha,
    sr=cs*inv_logit(100*cs),

    p=lambda+(1-2*lambda)*pnorm(beta*sr)
  ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='re')


pf_rel %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

pf_nm<-
  fit_nm$draws('mu',format = 'df') %>%
  mutate(
    alpha_1=exp(`mu[1]`),
    alpha_2=exp(`mu[2]`),
    alpha_3=exp(`mu[3]`),
    alpha_4=exp(`mu[4]`),
    alpha_5=exp(`mu[5]`),
    beta_1=exp(`mu[6]`),
    beta_2=exp(`mu[7]`),
    beta_3=exp(`mu[8]`),
    beta_4=exp(`mu[9]`),
    beta_5=exp(`mu[10]`),
    lambda=.5*inv_logit(`mu[11]`)
  ) %>% 
  pivot_longer(
    cols=contains(c('alpha','beta')),
    names_to = c("parameter", "at_idx"),
    names_sep = "_",
    values_to = "value"
    ) %>% 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) %>% 
  mutate(at_idx=as.double(at_idx)) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,at_idx=32-(30:34)+3,temperature=seq(28,35,.1),.draw=1:4000)) %>% 
  mutate(
    dfat=adapting-temperature,
    
    cs=dfat-alpha,
    sr=cs*inv_logit(100*cs),
    
    p=lambda+(1-2*lambda)*pnorm(beta*sr)
  ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='nm')

pf_nm %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

full_join(pf_abs,pf_rel) %>% 
  full_join(pf_nm) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=model,fill=model))+
  geom_line()+
  geom_ribbon(alpha=.5)+
  geom_vline(aes(xintercept = adapting),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  facet_grid(rows=vars(adapting))+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

#### Load fitted models ####
fit_abs <- readRDS("results/fits/1_2.rds")
fit_rel <- readRDS("results/fits/2_2.rds")
fit_nm <- readRDS("results/fits/3_2.rds")

#Extract posterior and reconstruct PF
pf_abs<-
  fit_abs$draws('mu',format = 'df') %>%
  mutate(
    rho=`mu[1]`,
    alpha=exp(`mu[2]`),
    beta=exp(`mu[3]`),
    lambda=.5*inv_logit(`mu[4]`)
  ) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,temperature=seq(29,36,.1),.draw=1:4000)) %>% 
  mutate(
    catt=temperature-baseline,
    caat=adapting-baseline,
    
    cs=catt-rho,
    sr=cs*inv_logit(100*cs),
    asr=sr*inv_logit(100*(sr-(caat+alpha-rho))),
    
    p=lambda+(1-2*lambda)*pnorm(beta*asr)
  ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='ab')

pf_abs %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

pf_rel<-
  fit_rel$draws('mu',format = 'df') %>%
  mutate(
    alpha=exp(`mu[1]`),
    beta=exp(`mu[2]`),
    lambda=.5*inv_logit(`mu[3]`)
  ) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,temperature=seq(29,36,.1),.draw=1:4000)) %>% 
  mutate(
    dfat=temperature-adapting,
    
    cs=dfat-alpha,
    sr=cs*inv_logit(100*cs),
    
    p=lambda+(1-2*lambda)*pnorm(beta*sr)
  ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='re')

pf_rel %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

pf_nm<-
  fit_nm$draws('mu',format = 'df') %>%
  mutate(
    alpha_1=exp(`mu[1]`),
    alpha_2=exp(`mu[2]`),
    alpha_3=exp(`mu[3]`),
    alpha_4=exp(`mu[4]`),
    alpha_5=exp(`mu[5]`),
    beta_1=exp(`mu[6]`),
    beta_2=exp(`mu[7]`),
    beta_3=exp(`mu[8]`),
    beta_4=exp(`mu[9]`),
    beta_5=exp(`mu[10]`),
    lambda=.5*inv_logit(`mu[11]`)
  ) %>% 
  pivot_longer(
    cols=contains(c('alpha','beta')),
    names_to = c("parameter", "at_idx"),
    names_sep = "_",
    values_to = "value"
  ) %>% 
  pivot_wider(
    names_from = parameter,
    values_from = value
  ) %>% 
  mutate(at_idx=as.double(at_idx)) %>% 
  full_join(expand.grid(baseline=32,adapting=30:34,at_idx=32-(30:34)+3,temperature=seq(29,36,.1),.draw=1:4000)) %>% 
  mutate(
    dfat=temperature-adapting,
    
    cs=dfat-alpha,
    sr=cs*inv_logit(100*cs),
    
    p=lambda+(1-2*lambda)*pnorm(beta*sr)
  ) %>% 
  group_by(temperature,adapting) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,.025),
    ub=quantile(p,.975)
  ) %>%
  mutate(model='nm')

pf_nm %>% 
  mutate(adapting=factor(adapting)) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=adapting,fill=adapting))+
  geom_line()+
  geom_ribbon(,alpha=.5)+
  geom_vline(aes(xintercept = 30),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 31),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 32),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 33),color='grey',linetype='dashed')+
  geom_vline(aes(xintercept = 34),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()

full_join(pf_abs,pf_rel) %>% 
  full_join(pf_nm) %>% 
  ggplot(aes(x=temperature,y=m,ymin=lb,ymax=ub,color=model,fill=model))+
  geom_line()+
  geom_ribbon(alpha=.5)+
  geom_vline(aes(xintercept = adapting),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = .5),color='grey',linetype='dashed')+
  geom_hline(aes(yintercept = 1),color='grey',linetype='dashed')+
  facet_grid(rows=vars(adapting))+
  labs(
    x='Temperature',
    y='P(correct response)',
    color='Adapting',
    fill='Adapting'
  )+
  theme_classic()
  
