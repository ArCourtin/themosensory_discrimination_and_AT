# Script used to plot the group-mean discrimination/rating psychometric functions of all four
# fitted models (1=absolute personal-baseline, 2=absolute fixed-reference, 3=relative,
# 4=non-mechanistic), one figure per domain x model, with cold/warm tasks overlaid (color=task)
# across the five adapting temperatures and the raw group-level data plotted underneath. Curve
# reconstruction math for each stan_models/*.stan file (ribbon/median style, nested-CI-via-alpha)
# follows plot_winning_models.R, extended here to loop over all four models instead of only the
# LOO winner.
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

n_draws<-1000
thin_draws<-function(df){
  if (nrow(df) > n_draws) slice_sample(df, n = n_draws) else df
}

#### Shared CI summary (nested 60/80/90/95% CI + median) ####
summarise_theta<-function(df){
  df %>%
    group_by(at,x) %>%
    summarise(
      lb_ci_95=quantile(theta,.025), ub_ci_95=quantile(theta,.975),
      lb_ci_90=quantile(theta,.05),  ub_ci_90=quantile(theta,.95),
      lb_ci_80=quantile(theta,.1),   ub_ci_80=quantile(theta,.9),
      lb_ci_60=quantile(theta,.2),   ub_ci_60=quantile(theta,.8),
      med_ci_50=quantile(theta,.5),
      .groups="drop"
    ) %>%
    pivot_longer(
      cols=matches("(lb|ub|med)_ci_\\d+"),
      names_to=c("stat","ci"),
      names_pattern="(lb|ub|med)_ci_(\\d+)",
      values_to="value"
    )
}

#### Group-mean curve reconstruction, one function per stan_models/*.stan file ####
gm_discrimination_absolute<-function(fit,at_seq,d){
  grid<-expand_grid(x=d,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,rho=`mu[1]`,beta=exp(`mu[2]`),lambda=.5*inv_logit(`mu[3]`),kappa=`mu[4]`) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      interval_sign=sign(x),
      target=at+abs(x),
      cx_target=target-nominal_baseline-rho,
      absolute_reading=cx_target*inv_logit(100*cx_target),
      mask_gate=inv_logit(100*(absolute_reading-(at-nominal_baseline-rho))),
      stim_rep=absolute_reading*mask_gate,
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    summarise_theta()
}

gm_discrimination_relative<-function(fit,at_seq,d){
  grid<-expand_grid(x=d,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,beta=exp(`mu[1]`),lambda=.5*inv_logit(`mu[2]`),kappa=`mu[3]`) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      interval_sign=sign(x),
      cx=abs(x),
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    summarise_theta()
}

gm_discrimination_non_mechanistic<-function(fit,at_seq,d){
  grid<-expand_grid(x=d,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(
      .draw,
      f_a3=`mu[1]`, f_a4=f_a3+`mu[2]`, f_a5=f_a4+`mu[3]`, f_a2=f_a3+`mu[4]`, f_a1=f_a2+`mu[5]`,
      f_b3=`mu[6]`, f_b4=f_b3+`mu[7]`, f_b5=f_b4+`mu[8]`, f_b2=f_b3+`mu[9]`, f_b1=f_b2+`mu[10]`,
      lambda=.5*inv_logit(`mu[11]`),
      kappa=`mu[12]`
    ) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      c=at-29,
      log_alpha=case_when(c==1~f_a1,c==2~f_a2,c==3~f_a3,c==4~f_a4,c==5~f_a5),
      log_beta =case_when(c==1~f_b1,c==2~f_b2,c==3~f_b3,c==4~f_b4,c==5~f_b5),
      interval_sign=sign(x),
      alpha=exp(log_alpha),
      beta=exp(log_beta),
      cx=abs(x)-alpha,
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa)
    ) %>%
    summarise_theta()
}

gm_rating_absolute<-function(fit,at_seq,x_grid,is_cold){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,intercept=`mu[1]`,slope=exp(`mu[2]`)) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(cx=if(is_cold) nominal_baseline-x else x-nominal_baseline,theta=inv_logit(intercept+slope*cx)) %>%
    summarise_theta()
}

gm_rating_relative<-function(fit,at_seq,x_grid,is_cold){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,intercept=`mu[1]`,slope=exp(`mu[2]`)) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(cx=if(is_cold) at-x else x-at,theta=inv_logit(intercept+slope*cx)) %>%
    summarise_theta()
}

gm_rating_non_mechanistic<-function(fit,at_seq,x_grid,is_cold){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(
      .draw,
      f_i3=`mu[1]`, f_i4=f_i3+`mu[2]`, f_i5=f_i4+`mu[3]`, f_i2=f_i3+`mu[4]`, f_i1=f_i2+`mu[5]`,
      f_s3=`mu[6]`, f_s4=f_s3+`mu[7]`, f_s5=f_s4+`mu[8]`, f_s2=f_s3+`mu[9]`, f_s1=f_s2+`mu[10]`
    ) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      c=at-29,
      intercept=case_when(c==1~f_i1,c==2~f_i2,c==3~f_i3,c==4~f_i4,c==5~f_i5),
      log_slope =case_when(c==1~f_s1,c==2~f_s2,c==3~f_s3,c==4~f_s4,c==5~f_s5),
      slope=exp(log_slope),
      cx=if(is_cold) at-x else x-at,
      theta=inv_logit(intercept+slope*cx)
    ) %>%
    summarise_theta()
}

# Models 1 (absolute, personal-baseline reference) and 2 (absolute, fixed reference) share the
# same stan model/parameterization - only the data fed to recorded_baseline_temperature differs -
# so the same group-mean reconstruction function applies to both.
gm_functions<-list(
  discrimination=list(gm_discrimination_absolute,gm_discrimination_absolute,gm_discrimination_relative,gm_discrimination_non_mechanistic),
  rating=list(gm_rating_absolute,gm_rating_absolute,gm_rating_relative,gm_rating_non_mechanistic)
)
model_labels<-c(`1`="absolute (personal baseline)",`2`="absolute (fixed reference)",`3`="relative",`4`="non-mechanistic")

task_labels<-c(`1`="cd",`2`="wd")
task_colors<-c(cd='#56B4E9',wd='#E69F00')

#### Raw data for overlay ####
# Discrimination: same filters/derivations as fit_discrimination_models_SLURM.R; x is the signed
# deviation of the recorded target temperature from the adapting temperature (matches the model's
# x-grid), binned per (task, AT condition, deviation) into a response proportion, as in
# visual_check_of_the_data.R's group-level summary. `at` uses the nominal baseline (32), matching
# the facet grid built from at_seq, not each participant's own recorded baseline.
raw_discrimination<-
  read_csv("data/d_at_2ifc_af.csv") %>%
  filter(baseline_flag==0, deviation_flag==0) %>%
  mutate(
    relative_adapting_temperature=round(adapting-baseline),
    at=nominal_baseline+relative_adapting_temperature,
    x=abs(recorded_temperature-adapting)*(active_interval-1.5)*2,
    chose_second=as.integer(if_else(active_interval==2,accuracy,1L-accuracy)),
    task_file_idx=task+1
  ) %>%
  group_by(task_file_idx,at,x) %>%
  summarise(m=mean(chose_second),n=sum(!is.na(chose_second)),.groups="drop") %>%
  mutate(task=task_labels[as.character(task_file_idx)])

# Rating: raw per-trial ratings (no binning, matching visual_check_of_the_data.R), x is the
# recorded absolute target temperature.
raw_rating<-
  read_csv("data/d_at_ratings_af.csv") %>%
  filter(baseline_flag==0, deviation_flag==0, confirmed==1, !is.na(recorded_temperature)) %>%
  mutate(
    relative_adapting_temperature=round(adapting-baseline),
    at=nominal_baseline+relative_adapting_temperature,
    x=recorded_temperature,
    rating=rating/100,
    task_file_idx=task+1,
    task=task_labels[as.character(task_file_idx)]
  )

#### Discrimination: all four models, both tasks overlaid (color = task) ####
d_grid<-list(`1`=d_cold,`2`=d_warm)

for(m in 1:4){
  curves<-map_dfr(names(task_labels),function(task){
    fit<-readRDS(paste0("results/fits/discrimination_",m,"_",task,".rds"))
    gm_functions$discrimination[[m]](fit,at_seq,d_grid[[task]]) %>%
      mutate(task=task_labels[[task]])
  })

  ribbons<-curves %>% filter(stat %in% c("lb","ub")) %>%
    pivot_wider(names_from="stat",values_from="value")
  medians<-curves %>% filter(stat=="med")

  p<-ribbons %>%
    ggplot()+
    geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task))+
    geom_line(data=medians,aes(x=x,y=value,color=task))+
    geom_point(data=raw_discrimination,aes(x=x,y=m,size=n,color=task),alpha=.4,shape=1)+
    geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
    geom_hline(aes(yintercept=.5),linetype='dotted',alpha=.5)+
    geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
    geom_vline(aes(xintercept=0),linetype='dashed',alpha=.5)+
    scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
    scale_fill_manual(values=task_colors,guide='none')+
    scale_color_manual(values=task_colors,guide='none')+
    scale_size_continuous(range=c(.5,3),guide='none')+
    scale_x_continuous(breaks=seq(-8,8,4))+
    scale_y_continuous(breaks=c(0,.5,1))+
    labs(
      alpha='',
      x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
      y='P(chose 2nd interval)',
      title=paste0("Discrimination - ",model_labels[[as.character(m)]])
    )+
    theme_minimal()+
    facet_grid(rows=vars(task),cols=vars(at))

  ggsave(paste0('figures/gm_discrimination_',m,'.png'),plot=p,units='cm',width=18,height=14)
}

#### Rating: all four models, both tasks overlaid (color = task) ####
task_is_cold<-c(`1`=TRUE,`2`=FALSE)
x_rating_grid<-list(`1`=x_rating_cold,`2`=x_rating_warm)

for(m in 1:4){
  curves<-map_dfr(names(task_labels),function(task){
    fit<-readRDS(paste0("results/fits/rating_",m,"_",task,".rds"))
    gm_functions$rating[[m]](fit,at_seq,x_rating_grid[[task]],task_is_cold[[task]]) %>%
      mutate(task=task_labels[[task]])
  })

  ribbons<-curves %>% filter(stat %in% c("lb","ub")) %>%
    pivot_wider(names_from="stat",values_from="value")
  medians<-curves %>% filter(stat=="med")

  p<-ribbons %>%
    ggplot()+
    geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task))+
    geom_line(data=medians,aes(x=x,y=value,color=task))+
    geom_point(data=raw_rating,aes(x=x,y=rating,color=task),alpha=.1,shape=16,size=.6)+
    geom_hline(aes(yintercept=0),linetype='dotted',alpha=.5)+
    geom_hline(aes(yintercept=1),linetype='dotted',alpha=.5)+
    scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
    scale_fill_manual(values=task_colors,guide='none')+
    scale_color_manual(values=task_colors,guide='none')+
    scale_x_continuous(breaks=seq(22,45,5))+
    scale_y_continuous(breaks=c(0,.5,1))+
    labs(
      alpha='',
      x='Absolute target temperature',
      y='Rating',
      title=paste0("Rating - ",model_labels[[as.character(m)]])
    )+
    theme_minimal()+
    facet_grid(rows=vars(task),cols=vars(at))

  ggsave(paste0('figures/gm_rating_',m,'.png'),plot=p,units='cm',width=18,height=14)
}
