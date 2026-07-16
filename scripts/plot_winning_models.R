# Script used to plot the group-mean discrimination/rating functions of the LOO-winning model,
# separately for the cold and warm tasks, across the five adapting temperatures (baseline 32,
# relative AT in -2:2). Ribbon/median style and nested-CI-via-alpha follow plot_priors.R; unlike
# plot_priors.R (which facets subject- vs group-level draws), only the group-mean curve is shown
# here, so that facet row is repurposed to color-code the two tasks in the same panels instead.
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(loo)

rm(list=ls())

inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

baseline<-32
at_seq<-baseline+(-2:2)
d<-seq(-13,13,.1)          # discrimination: signed target deviation
x_rating<-seq(30,45,.1)    # rating: absolute target temperature

n_draws<-1000
thin_draws<-function(df){
  if (nrow(df) > n_draws) slice_sample(df, n = n_draws) else df
}

#### Shared CI summary (nested 60/80/90/95% CI + median, matching plot_priors.R) ####
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
      cx_target=target-baseline-rho,
      absolute_reading=cx_target*inv_logit(100*cx_target),
      mask_gate=inv_logit(100*(absolute_reading-(at-baseline-rho))),
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

gm_rating_absolute<-function(fit,at_seq,x_grid){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,intercept=`mu[1]`,slope=exp(`mu[2]`)) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(cx=x-baseline,theta=inv_logit(intercept+slope*cx)) %>%
    summarise_theta()
}

gm_rating_relative<-function(fit,at_seq,x_grid){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,intercept=`mu[1]`,slope=exp(`mu[2]`)) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(cx=x-at,theta=inv_logit(intercept+slope*cx)) %>%
    summarise_theta()
}

gm_rating_non_mechanistic<-function(fit,at_seq,x_grid){
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
      cx=x-at,
      theta=inv_logit(intercept+slope*cx)
    ) %>%
    summarise_theta()
}

gm_functions<-list(
  discrimination=list(gm_discrimination_absolute,gm_discrimination_relative,gm_discrimination_non_mechanistic),
  rating=list(gm_rating_absolute,gm_rating_relative,gm_rating_non_mechanistic)
)

#### Determine the LOO-winning model per domain and task ####
# results/loo/{domain}_{model}_{task}.rds: model 1=absolute, 2=relative, 3=non-mechanistic;
# task 1=cold, task 2=warm. results/fits/ mirrors the same naming.
winning_model_idx<-function(domain,task){
  loos<-list(
    absolute=readRDS(paste0("results/loo/",domain,"_1_",task,".rds")),
    relative=readRDS(paste0("results/loo/",domain,"_2_",task,".rds")),
    `non-mechanistic`=readRDS(paste0("results/loo/",domain,"_3_",task,".rds"))
  )
  winner<-rownames(loo_compare(loos))[1]
  match(winner,c("absolute","relative","non-mechanistic"))
}

task_labels<-c(`1`="cd",`2`="wd")
# Fixed per-task color, matching the palette used project-wide in AUD_thresholding.
task_colors<-c(cd='#56B4E9',wd='#E69F00')

#### Discrimination: winning model per task, both tasks overlaid (color = task) ####
discrimination_curves<-map_dfr(names(task_labels),function(task){
  m<-winning_model_idx("discrimination",task)
  fit<-readRDS(paste0("results/fits/discrimination_",m,"_",task,".rds"))
  gm_functions$discrimination[[m]](fit,at_seq,d) %>%
    mutate(task=task_labels[[task]])
})

ribbons<-discrimination_curves %>% filter(stat %in% c("lb","ub")) %>%
  pivot_wider(names_from="stat",values_from="value")
medians<-discrimination_curves %>% filter(stat=="med")

ribbons %>%
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task)
  )+
  geom_line(
    data=medians,
    aes(x=x,y=value,color=task)
  )+
  geom_hline(
    aes(yintercept=.5),
    linetype='dotted',
    alpha=.5
  )+
  geom_hline(
    aes(yintercept=1),
    linetype='dotted',
    alpha=.5
  )+
  geom_vline(aes(xintercept=0),linetype='dotted',alpha=.5)+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_fill_manual(values=task_colors,guide='none')+
  scale_color_manual(values=task_colors,guide='none')+
  scale_x_continuous(breaks=seq(-12,12,4))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
    y='P(chose 2nd interval)'
  )+
  theme_minimal()+
  facet_grid(rows=vars(task),cols=vars(at))

ggsave('figures/gm_discrimination.png',units='cm',width=18,height=14)

#### Rating: winning model per task, both tasks overlaid (color = task) ####
rating_curves<-map_dfr(names(task_labels),function(task){
  m<-winning_model_idx("rating",task)
  fit<-readRDS(paste0("results/fits/rating_",m,"_",task,".rds"))
  gm_functions$rating[[m]](fit,at_seq,x_rating) %>%
    mutate(task=task_labels[[task]])
})

ribbons<-rating_curves %>% filter(stat %in% c("lb","ub")) %>%
  pivot_wider(names_from="stat",values_from="value")
medians<-rating_curves %>% filter(stat=="med")

ribbons %>%
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task)
  )+
  geom_line(
    data=medians,
    aes(x=x,y=value,color=task)
  )+
  geom_hline(
    aes(yintercept=0),
    linetype='dotted',
    alpha=.5
  )+
  geom_hline(
    aes(yintercept=1),
    linetype='dotted',
    alpha=.5
  )+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_fill_manual(values=task_colors,guide='none')+
  scale_color_manual(values=task_colors,guide='none')+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='Mean rating'
  )+
  theme_minimal()+
  facet_grid(rows=vars(task),cols=vars(at))

ggsave('figures/gm_rating.png',units='cm',width=18,height=14)
