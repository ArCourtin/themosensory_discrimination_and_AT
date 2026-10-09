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
d_cold<-seq(-2,2,.1)               # discrimination: signed target deviation, cold task
d_warm<-seq(-8,8,.1)               # discrimination: signed target deviation, warm task
x_rating_cold<-seq(22,32,.1)       # rating: absolute target temperature, cold task
x_rating_warm<-seq(32,45,.1)       # rating: absolute target temperature, warm task

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
    transmute(
      .draw,xi=exp(`mu[1]`),beta=exp(`mu[2]`),
      f_l3=`mu[3]`, f_l4=f_l3+`mu[4]`, f_l5=f_l4+`mu[5]`, f_l2=f_l3+`mu[6]`, f_l1=f_l2+`mu[7]`,
      f_k3=`mu[8]`, f_k4=f_k3+`mu[9]`, f_k5=f_k4+`mu[10]`, f_k2=f_k3+`mu[11]`, f_k1=f_k2+`mu[12]`
    ) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      c=at-29,
      lambda=.5*inv_logit(pick_condition(c,f_l1,f_l2,f_l3,f_l4,f_l5)),
      kappa=pick_condition(c,f_k1,f_k2,f_k3,f_k4,f_k5),
      interval_sign=sign(x),
      target=at+abs(x),
      cx_target=target-baseline-xi,
      absolute_reading=cx_target*inv_logit(100*cx_target),
      mask_gate=inv_logit(100*(absolute_reading-(at-baseline-xi))),
      stim_rep=absolute_reading*mask_gate,
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa)
    ) %>%
    summarise_theta()
}

gm_discrimination_relative<-function(fit,at_seq,d){
  grid<-expand_grid(x=d,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(
      .draw,xi=exp(`mu[1]`),beta=exp(`mu[2]`),
      f_l3=`mu[3]`, f_l4=f_l3+`mu[4]`, f_l5=f_l4+`mu[5]`, f_l2=f_l3+`mu[6]`, f_l1=f_l2+`mu[7]`,
      f_k3=`mu[8]`, f_k4=f_k3+`mu[9]`, f_k5=f_k4+`mu[10]`, f_k2=f_k3+`mu[11]`, f_k1=f_k2+`mu[12]`
    ) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      c=at-29,
      lambda=.5*inv_logit(pick_condition(c,f_l1,f_l2,f_l3,f_l4,f_l5)),
      kappa=pick_condition(c,f_k1,f_k2,f_k3,f_k4,f_k5),
      interval_sign=sign(x),
      cx=abs(x)-xi,
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa)
    ) %>%
    summarise_theta()
}

gm_discrimination_non_mechanistic<-function(fit,at_seq,d){
  grid<-expand_grid(x=d,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(
      .draw,
      f_x3=`mu[1]`, f_x4=f_x3+`mu[2]`, f_x5=f_x4+`mu[3]`, f_x2=f_x3+`mu[4]`, f_x1=f_x2+`mu[5]`,
      f_b3=`mu[6]`, f_b4=f_b3+`mu[7]`, f_b5=f_b4+`mu[8]`, f_b2=f_b3+`mu[9]`, f_b1=f_b2+`mu[10]`,
      f_l3=`mu[11]`, f_l4=f_l3+`mu[12]`, f_l5=f_l4+`mu[13]`, f_l2=f_l3+`mu[14]`, f_l1=f_l2+`mu[15]`,
      f_k3=`mu[16]`, f_k4=f_k3+`mu[17]`, f_k5=f_k4+`mu[18]`, f_k2=f_k3+`mu[19]`, f_k1=f_k2+`mu[20]`
    ) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(
      c=at-29,
      log_xi=case_when(c==1~f_x1,c==2~f_x2,c==3~f_x3,c==4~f_x4,c==5~f_x5),
      log_beta =case_when(c==1~f_b1,c==2~f_b2,c==3~f_b3,c==4~f_b4,c==5~f_b5),
      lambda=.5*inv_logit(pick_condition(c,f_l1,f_l2,f_l3,f_l4,f_l5)),
      kappa=pick_condition(c,f_k1,f_k2,f_k3,f_k4,f_k5),
      interval_sign=sign(x),
      xi=exp(log_xi),
      beta=exp(log_beta),
      cx=abs(x)-xi,
      stim_rep=cx*inv_logit(100*cx),
      theta=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa)
    ) %>%
    summarise_theta()
}

gm_rating_absolute<-function(fit,at_seq,x_grid,is_cold){
  grid<-expand_grid(x=x_grid,at=at_seq)
  fit$draws('mu',format='df') %>%
    transmute(.draw,intercept=`mu[1]`,slope=exp(`mu[2]`)) %>%
    thin_draws() %>%
    cross_join(grid) %>%
    mutate(cx=if(is_cold) baseline-x else x-baseline,theta=inv_logit(intercept+slope*cx)) %>%
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

gm_functions<-list(
  discrimination=list(`1`=gm_discrimination_absolute,`2`=gm_discrimination_relative,`3`=gm_discrimination_non_mechanistic),
  rating=list(`1`=gm_rating_absolute,`2`=gm_rating_relative,`3`=gm_rating_non_mechanistic)
)

#### Determine the LOO-winning model per domain and task ####
# results/loo/{domain}_{model}_{task}.rds: model 1=absolute, 2=relative, 3=non-mechanistic;
# task 1=cold, task 2=warm. results/fits/ mirrors the same naming.
model_labels<-c(`1`="absolute",`2`="relative",`3`="non-mechanistic")
winning_model_idx<-function(domain,task){
  loos<-map(names(model_labels),~readRDS(paste0("results/loo/",domain,"_",.x,"_",task,".rds")))
  names(loos)<-model_labels
  winner<-rownames(loo_compare(loos))[1]
  as.integer(names(model_labels)[match(winner,model_labels)])
}

task_labels<-c(`1`="cd",`2`="wd")
# Fixed per-task color, matching the palette used project-wide in AUD_thresholding.
task_colors<-c(cd='#56B4E9',wd='#E69F00')

#### Discrimination: winning model per task, both tasks overlaid (color = task) ####
d_grid<-list(`1`=d_cold,`2`=d_warm)

discrimination_curves<-map_dfr(names(task_labels),function(task){
  m<-winning_model_idx("discrimination",task)
  fit<-readRDS(paste0("results/fits/discrimination_",m,"_",task,".rds"))
  gm_functions$discrimination[[as.character(m)]](fit,at_seq,d_grid[[task]]) %>%
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
    aes(yintercept=0),
    linetype='dotted',
    alpha=.5
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
  geom_vline(aes(xintercept=0),linetype='dashed',alpha=.5)+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_fill_manual(values=task_colors,guide='none')+
  scale_color_manual(values=task_colors,guide='none')+
  scale_x_continuous(breaks=seq(-8,8,4))+
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
task_is_cold<-c(`1`=TRUE,`2`=FALSE)
x_rating_grid<-list(`1`=x_rating_cold,`2`=x_rating_warm)

rating_curves<-map_dfr(names(task_labels),function(task){
  m<-winning_model_idx("rating",task)
  fit<-readRDS(paste0("results/fits/rating_",m,"_",task,".rds"))
  gm_functions$rating[[as.character(m)]](fit,at_seq,x_rating_grid[[task]],task_is_cold[[task]]) %>%
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
  scale_x_continuous(breaks=seq(22,45,5))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='Mean rating'
  )+
  theme_minimal()+
  facet_grid(rows=vars(task),cols=vars(at))

ggsave('figures/gm_rating.png',units='cm',width=18,height=14)
