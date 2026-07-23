# Script used to plot each participant's individual discrimination/rating psychometric function
# against their own raw data, one figure per domain x model x task (cold/warm), faceted by
# participant (rows) and adapting-temperature condition (columns). Unlike make_gm_plots.R (which
# reconstructs the group-mean curve with credible ribbons from posterior `mu` draws), this script
# reconstructs one point-estimate curve per participant from the posterior median of each
# participant's own transformed parameters (e.g. `rho[p]`, `beta[p]`) - no ribbons, since raw
# per-participant data is too sparse for a per-draw uncertainty band to be informative here.
# Curve reconstruction math for each stan_models/*.stan file otherwise follows the same formulas
# as make_gm_plots.R / plot_winning_models.R, indexed by participant instead of by draw. The
# absolute model's individual curves use each participant's own recorded pre-session baseline
# temperature (model 1) or the fixed common reference of 32 (model 2), matching what each model
# was actually fit with, rather than the nominal baseline used for the group-mean grid.
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Written with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)

rm(list=ls())

inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

nominal_baseline<-32
at_seq<-nominal_baseline+(-2:2)
d_cold<-seq(-2,2,.1)               # discrimination: signed target deviation, cold task
d_warm<-seq(-8,8,.1)               # discrimination: signed target deviation, warm task
x_rating_cold<-seq(22,32,.1)       # rating: absolute target temperature, cold task
x_rating_warm<-seq(32,45,.1)       # rating: absolute target temperature, warm task

model_labels<-c(`1`="absolute (personal baseline)",`2`="absolute (fixed reference)",`3`="relative",`4`="non-mechanistic")
task_labels<-c(`1`="cd",`2`="wd")

#### Per-participant parameter extraction from posterior medians ####
# Vector[P] transformed parameters (e.g. rho[p], beta[p]) - one row per participant per variable.
ind_params_vecP<-function(fit,vars){
  fit$summary(vars) %>%
    transmute(
      base=str_extract(variable,"^[^\\[]+"),
      participant=as.integer(str_extract(variable,"(?<=\\[)\\d+(?=\\])")),
      median
    ) %>%
    pivot_wider(names_from=base,values_from=median)
}

# Matrix[C,P] transformed parameters (e.g. alpha[c,p]) - one row per (condition, participant) per variable.
ind_params_matCP<-function(fit,vars){
  fit$summary(vars) %>%
    mutate(
      base=str_extract(variable,"^[^\\[]+"),
      idx=str_extract(variable,"(?<=\\[)[^\\]]+(?=\\])")
    ) %>%
    separate(idx,into=c("c","participant"),sep=",",convert=TRUE) %>%
    select(base,c,participant,median) %>%
    pivot_wider(names_from=base,values_from=median)
}

#### Individual-curve reconstruction, one function per stan_models/*.stan file ####
# participant_baseline (only used by the absolute model) is a df(participant,baseline_i) giving
# each participant's own reference temperature - their recorded pre-session baseline for model 1,
# a fixed 32 for model 2. The unused argument is kept in the other models' signatures so all four
# can be called uniformly from the same loop.
ind_discrimination_absolute<-function(fit,at_seq,d,participant_baseline){
  grid<-expand_grid(x=d,at=at_seq)
  ind_params_vecP(fit,c('rho','beta','lambda','kappa')) %>%
    left_join(participant_baseline,by='participant') %>%
    cross_join(grid) %>%
    mutate(
      interval_sign=sign(x),
      target=at+abs(x),
      cx_target=target-baseline_i-rho,
      absolute_reading=cx_target*inv_logit(100*cx_target),
      mask_gate=inv_logit(100*(absolute_reading-(at-baseline_i-rho))),
      stim_rep=absolute_reading*mask_gate,
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    select(participant,at,x,theta)
}

ind_discrimination_relative<-function(fit,at_seq,d,participant_baseline=NULL){
  grid<-expand_grid(x=d,at=at_seq)
  ind_params_vecP(fit,c('beta','lambda','kappa')) %>%
    cross_join(grid) %>%
    mutate(
      interval_sign=sign(x),
      cx=abs(x),
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    select(participant,at,x,theta)
}

ind_discrimination_non_mechanistic<-function(fit,at_seq,d,participant_baseline=NULL){
  grid<-expand_grid(x=d,at=at_seq) %>% mutate(c=at-29)
  ab<-ind_params_matCP(fit,c('alpha','beta'))
  lk<-ind_params_vecP(fit,c('lambda','kappa'))
  ab %>%
    inner_join(grid,by='c') %>%
    left_join(lk,by='participant') %>%
    mutate(
      interval_sign=sign(x),
      cx=abs(x)-alpha,
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    select(participant,at,x,theta)
}

ind_rating_absolute<-function(fit,at_seq,x_grid,is_cold,participant_baseline){
  grid<-expand_grid(x=x_grid,at=at_seq)
  ind_params_vecP(fit,c('intercept','slope')) %>%
    left_join(participant_baseline,by='participant') %>%
    cross_join(grid) %>%
    mutate(cx=if(is_cold) baseline_i-x else x-baseline_i,theta=inv_logit(intercept+slope*cx)) %>%
    select(participant,at,x,theta)
}

ind_rating_relative<-function(fit,at_seq,x_grid,is_cold,participant_baseline=NULL){
  grid<-expand_grid(x=x_grid,at=at_seq)
  ind_params_vecP(fit,c('intercept','slope')) %>%
    cross_join(grid) %>%
    mutate(cx=if(is_cold) at-x else x-at,theta=inv_logit(intercept+slope*cx)) %>%
    select(participant,at,x,theta)
}

ind_rating_non_mechanistic<-function(fit,at_seq,x_grid,is_cold,participant_baseline=NULL){
  grid<-expand_grid(x=x_grid,at=at_seq) %>% mutate(c=at-29)
  ind_params_matCP(fit,c('intercept','slope')) %>%
    inner_join(grid,by='c') %>%
    mutate(cx=if(is_cold) at-x else x-at,theta=inv_logit(intercept+slope*cx)) %>%
    select(participant,at,x,theta)
}

ind_functions<-list(
  discrimination=list(ind_discrimination_absolute,ind_discrimination_absolute,ind_discrimination_relative,ind_discrimination_non_mechanistic),
  rating=list(ind_rating_absolute,ind_rating_absolute,ind_rating_relative,ind_rating_non_mechanistic)
)

fig_height<-function(P){ max(10,2.2*P) }

#### Discrimination ####
# Same filters/derivations/participant-renumbering as fit_discrimination_models_SLURM.R.
discrimination_data<-
  read_csv("data/d_at_2ifc_af.csv") %>%
  filter(baseline_flag==0, deviation_flag==0) %>%
  mutate(
    relative_adapting_temperature=round(adapting-baseline),
    adapting_temperature_idx=3+relative_adapting_temperature,
    interval_sign=if_else(active_interval==2,1,-1),
    chose_second=as.integer(if_else(active_interval==2,accuracy,1L-accuracy))
  )
participant<-unique(discrimination_data$participant)
for(pdx in seq_along(participant)){
  discrimination_data$participant[discrimination_data$participant==participant[pdx]]<-pdx
}

d_grid<-list(`1`=d_cold,`2`=d_warm)

for(task in names(task_labels)){
  t<-as.integer(task)-1
  sample_data<-discrimination_data %>% filter(task==t)
  P<-n_distinct(sample_data$participant)

  personal_baseline_df<-sample_data %>% distinct(participant,baseline) %>% rename(baseline_i=baseline)
  fixed_baseline_df<-personal_baseline_df %>% mutate(baseline_i=nominal_baseline)
  baseline_by_model<-list(`1`=personal_baseline_df,`2`=fixed_baseline_df,`3`=NULL,`4`=NULL)

  raw_ind<-sample_data %>%
    mutate(
      at=nominal_baseline+relative_adapting_temperature,
      x=abs(recorded_temperature-adapting)*(active_interval-1.5)*2
    ) %>%
    group_by(participant,at,x) %>%
    summarise(m=mean(chose_second),n=sum(!is.na(chose_second)),.groups="drop")

  for(m in 1:4){
    fit<-readRDS(paste0("results/fits/discrimination_",m,"_",task,".rds"))
    curves<-ind_functions$discrimination[[m]](fit,at_seq,d_grid[[task]],baseline_by_model[[as.character(m)]])

    p<-curves %>%
      ggplot()+
      geom_line(aes(x=x,y=theta))+
      geom_point(data=raw_ind,aes(x=x,y=m,alpha=n),size=.8)+
      geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
      geom_hline(aes(yintercept=.5),linetype='dotted',alpha=.5)+
      geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
      geom_vline(aes(xintercept=0),linetype='dashed',alpha=.5)+
      scale_alpha_continuous(range=c(.3,1),guide='none')+
      scale_x_continuous(breaks=seq(-8,8,4))+
      scale_y_continuous(breaks=c(0,.5,1))+
      labs(
        x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
        y='P(chose 2nd interval)',
        title=paste0("Discrimination (",task_labels[[task]],") - ",model_labels[[as.character(m)]])
      )+
      theme_minimal()+
      facet_grid(rows=vars(participant),cols=vars(at))

    ggsave(paste0('figures/individual_discrimination_',m,'_',task,'.png'),plot=p,units='cm',width=18,height=fig_height(P),limitsize=FALSE)
  }
}

#### Rating ####
# Same filters/derivations/participant-renumbering as fit_rating_models_SLURM.R.
rating_data<-
  read_csv("data/d_at_ratings_af.csv") %>%
  filter(baseline_flag==0, deviation_flag==0, confirmed==1, !is.na(recorded_temperature)) %>%
  mutate(
    relative_adapting_temperature=round(adapting-baseline),
    adapting_temperature_idx=3+relative_adapting_temperature
  )
participant<-unique(rating_data$participant)
for(pdx in seq_along(participant)){
  rating_data$participant[rating_data$participant==participant[pdx]]<-pdx
}

task_is_cold<-c(`1`=TRUE,`2`=FALSE)
x_rating_grid<-list(`1`=x_rating_cold,`2`=x_rating_warm)

for(task in names(task_labels)){
  t<-as.integer(task)-1
  sample_data<-rating_data %>% filter(task==t)
  P<-n_distinct(sample_data$participant)

  personal_baseline_df<-sample_data %>% distinct(participant,baseline) %>% rename(baseline_i=baseline)
  fixed_baseline_df<-personal_baseline_df %>% mutate(baseline_i=nominal_baseline)
  baseline_by_model<-list(`1`=personal_baseline_df,`2`=fixed_baseline_df,`3`=NULL,`4`=NULL)

  raw_ind<-sample_data %>%
    mutate(
      at=nominal_baseline+relative_adapting_temperature,
      x=recorded_temperature,
      rating=rating/100
    )

  for(m in 1:4){
    fit<-readRDS(paste0("results/fits/rating_",m,"_",task,".rds"))
    curves<-ind_functions$rating[[m]](fit,at_seq,x_rating_grid[[task]],task_is_cold[[task]],baseline_by_model[[as.character(m)]])

    p<-curves %>%
      ggplot()+
      geom_line(aes(x=x,y=theta))+
      geom_point(data=raw_ind,aes(x=x,y=rating),alpha=.3,size=.8)+
      geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
      geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
      scale_x_continuous(breaks=seq(22,45,5))+
      scale_y_continuous(breaks=c(0,.5,1))+
      labs(
        x='Absolute target temperature',
        y='Rating',
        title=paste0("Rating (",task_labels[[task]],") - ",model_labels[[as.character(m)]])
      )+
      theme_minimal()+
      facet_grid(rows=vars(participant),cols=vars(at))

    ggsave(paste0('figures/individual_rating_',m,'_',task,'.png'),plot=p,units='cm',width=18,height=fig_height(P),limitsize=FALSE)
  }
}
