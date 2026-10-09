#This script is used to illustrate the priors used for the different models
#Author: A.S. Courtin
#Licence: MIT
#Edited with the assistance of Claude Code (Anthropic).

library(tidyverse)

inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

M=10^3
idx=1:M
at=30:34
x=seq(30,45,.1)
grid_r<-expand.grid(x=x,idx=idx,at=at)
# Response-space axis for the discrimination models: x is the SIGNED deviation from the
# adapting temperature. |x| is the physical target deviation (target = at + |x|); the sign
# encodes which interval held the deviating stimulus (+ = 2nd interval, so interval_sign=+1).
d=seq(-13,13,.1)
grid_d<-expand.grid(x=d,idx=idx,at=at)

# Task conventions: cold (cd) blue, warm (wd) orange. Discrimination beta priors differ by task
# (matches the beta_mean/beta_sd in stan_models/discrimination_*.stan).
task_colors<-c(cd='#56B4E9',wd='#E69F00')
task_beta_prior<-list(
  cd=c(mean=0.552-30*0.0157-log(sqrt(2)),sd=sqrt(1.1^2+log(sqrt(2))^2)),
  wd=c(mean=-0.786-30*0.00239-log(sqrt(2)),sd=sqrt(0.865^2+log(sqrt(2))^2))
)

#### Discrimination - absolute coding####
for(task in names(task_colors)){
  beta_mean<-task_beta_prior[[task]][['mean']]
  beta_sd<-task_beta_prior[[task]][['sd']]

  discrimination_absolute<-
    tibble(idx=idx) %>% 
    rowwise() %>% 
    mutate(
      mu_log_xi=rnorm(1,-2,0.5),
      mu_log_beta=rnorm(1,beta_mean,beta_sd),
      mu_logit_lambda=rnorm(1,-3.85+30*0.00171,0.72),
      tau_log_xi=abs(rnorm(1,0,0.5)),
      tau_log_beta=abs(rnorm(1,0,beta_sd)),
      tau_logit_lambda=abs(rnorm(1,0,0.72)),
      xi=exp(rnorm(1,mu_log_xi,tau_log_xi)),
      beta=exp(rnorm(1,mu_log_beta,tau_log_beta)),
      lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
      mu_kappa=rnorm(1,0,0.5),
      tau_kappa=abs(rnorm(1,0,0.5)),
      kappa=rnorm(1,mu_kappa,tau_kappa)
    ) %>%
    ungroup() %>%
    full_join(grid_d) %>%
    mutate(
      # Evidence is zero at and below the adapting temperature and equal to the raw absolute
      # (baseline+xi anchored) reading above it - no separate detection threshold (matches
      # discrimination_absolute_coding.stan).
      interval_sign=sign(x),
      target=at+abs(x),
      cx_target_s=target-32-xi,
      absolute_reading_s=cx_target_s*inv_logit(cx_target_s*100),
      mask_gate_s=inv_logit(100*(absolute_reading_s-(at-32-xi))),
      stim_rep=absolute_reading_s*mask_gate_s,
      theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa),

      cx_target_g=target-32-exp(mu_log_xi),
      absolute_reading_g=cx_target_g*inv_logit(cx_target_g*100),
      mask_gate_g=inv_logit(100*(absolute_reading_g-(at-32-exp(mu_log_xi)))),
      stim_rep_g=absolute_reading_g*mask_gate_g,
      theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*exp(mu_log_beta)*stim_rep_g+mu_kappa)
      ) %>%
    group_by(at,x) %>% 
    summarise(
      s_lb_ci_95=quantile(theta_s,0.025),
      s_ub_ci_95=quantile(theta_s,0.975),
      s_lb_ci_90=quantile(theta_s,0.05),
      s_ub_ci_90=quantile(theta_s,0.95),
      s_lb_ci_80=quantile(theta_s,0.1),
      s_ub_ci_80=quantile(theta_s,0.9),
      s_lb_ci_60=quantile(theta_s,0.20),
      s_ub_ci_60=quantile(theta_s,0.80),
      g_lb_ci_95=quantile(theta_g,0.025),
      g_ub_ci_95=quantile(theta_g,0.975),
      g_lb_ci_90=quantile(theta_g,0.05),
      g_ub_ci_90=quantile(theta_g,0.95),
      g_lb_ci_80=quantile(theta_g,0.1),
      g_ub_ci_80=quantile(theta_g,0.9),
      g_lb_ci_60=quantile(theta_g,0.20),
      g_ub_ci_60=quantile(theta_g,0.80),
      s_med_ci_50 = quantile(theta_s, 0.5),
      g_med_ci_50 = quantile(theta_g, 0.5)
    ) %>% 
    pivot_longer(
      cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
      names_to = c("level","stat","ci"),
      names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
      values_to = "value"
    ) %>% 
    filter(TRUE) # response space: keep the full signed deviation range

  ribbons <- discrimination_absolute %>% filter(stat %in% c("lb","ub"))%>% 
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    )
  medians <- discrimination_absolute %>% filter(stat == "med")

  p <- ribbons %>%
    ggplot()+
    geom_ribbon(
      aes(x=x,ymin = lb,ymax=ub,alpha=ci),
      fill=task_colors[[task]]
    )+
    geom_line(
      data=medians,
      aes(x=x,y=value),
      color=task_colors[[task]]
    )+
    geom_hline(
      aes(yintercept = .5),
      linetype='dotted',
      alpha=.5
    )+
    geom_hline(
      aes(yintercept = 1),
      linetype='dotted',
      alpha=.5
    )+
    geom_vline(aes(xintercept=0),linetype='dotted',alpha=.5)+
    scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
    scale_x_continuous(breaks=seq(-12,12,4))+
    scale_y_continuous(breaks=c(0,.5,1))+
    labs(
      alpha='',
      x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
      y='P(chose 2nd interval)'
    )+
    theme_minimal()+
    facet_grid(rows=vars(level),cols=vars(at))

  ggsave(paste0('figures/d_a_ppp_',task,'.png'),plot=p,units = 'cm',width = 18,height = 10)
}

#### Discrimination - relative coding####
for(task in names(task_colors)){
  beta_mean<-task_beta_prior[[task]][['mean']]
  beta_sd<-task_beta_prior[[task]][['sd']]

  discrimination_relative<-
    tibble(idx=idx) %>% 
    rowwise() %>% 
    mutate(
      mu_log_xi=rnorm(1,-2,0.5),
      mu_log_beta=rnorm(1,beta_mean,beta_sd),
      mu_logit_lambda=rnorm(1,-3.85+30*0.00171,0.72),
      tau_log_xi=abs(rnorm(1,0,0.5)),
      tau_log_beta=abs(rnorm(1,0,beta_sd)),
      tau_logit_lambda=abs(rnorm(1,0,0.72)),
      xi=exp(rnorm(1,mu_log_xi,tau_log_xi)),
      beta=exp(rnorm(1,mu_log_beta,tau_log_beta)),
      lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
      mu_kappa=rnorm(1,0,0.5),
      tau_kappa=abs(rnorm(1,0,0.5)),
      kappa=rnorm(1,mu_kappa,tau_kappa)
    ) %>%
    ungroup() %>%
    full_join(grid_d) %>%
    mutate(
      # Evidence is the soft-rectified deviation from the adapting temperature beyond the
      # hard threshold xi, i.e. the silent range (matches discrimination_relative_coding.stan).
      interval_sign=sign(x),
      cx_s=abs(x)-xi,
      stim_rep=cx_s*inv_logit(cx_s*100),
      theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa),

      cx_g=abs(x)-exp(mu_log_xi),
      stim_rep_g=cx_g*inv_logit(cx_g*100),
      theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*exp(mu_log_beta)*stim_rep_g+mu_kappa)
    ) %>%
    group_by(at,x) %>% 
    summarise(
      s_lb_ci_95=quantile(theta_s,0.025),
      s_ub_ci_95=quantile(theta_s,0.975),
      s_lb_ci_90=quantile(theta_s,0.05),
      s_ub_ci_90=quantile(theta_s,0.95),
      s_lb_ci_80=quantile(theta_s,0.1),
      s_ub_ci_80=quantile(theta_s,0.9),
      s_lb_ci_60=quantile(theta_s,0.20),
      s_ub_ci_60=quantile(theta_s,0.80),
      g_lb_ci_95=quantile(theta_g,0.025),
      g_ub_ci_95=quantile(theta_g,0.975),
      g_lb_ci_90=quantile(theta_g,0.05),
      g_ub_ci_90=quantile(theta_g,0.95),
      g_lb_ci_80=quantile(theta_g,0.1),
      g_ub_ci_80=quantile(theta_g,0.9),
      g_lb_ci_60=quantile(theta_g,0.20),
      g_ub_ci_60=quantile(theta_g,0.80),
      s_med_ci_50 = quantile(theta_s, 0.5),
      g_med_ci_50 = quantile(theta_g, 0.5)
    ) %>% 
    pivot_longer(
      cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
      names_to = c("level","stat","ci"),
      names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
      values_to = "value"
    )%>% 
    filter(TRUE) # response space: keep the full signed deviation range

  ribbons <- discrimination_relative %>% filter(stat %in% c("lb","ub"))%>% 
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    )
  medians <- discrimination_relative %>% filter(stat == "med")

  p <- ribbons %>%
    ggplot()+
    geom_ribbon(
      aes(x=x,ymin = lb,ymax=ub,alpha=ci),
      fill=task_colors[[task]]
    )+
    geom_line(
      data=medians,
      aes(x=x,y=value),
      color=task_colors[[task]]
    )+
    geom_hline(
      aes(yintercept = .5),
      linetype='dotted',
      alpha=.5
    )+
    geom_hline(
      aes(yintercept = 1),
      linetype='dotted',
      alpha=.5
    )+
    geom_vline(aes(xintercept=0),linetype='dotted',alpha=.5)+
    scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
    scale_x_continuous(breaks=seq(-12,12,4))+
    scale_y_continuous(breaks=c(0,.5,1))+
    labs(
      alpha='',
      x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
      y='P(chose 2nd interval)'
    )+
    theme_minimal()+
    facet_grid(rows=vars(level),cols=vars(at))

  ggsave(paste0('figures/d_r_ppp_',task,'.png'),plot=p,units = 'cm',width = 18,height = 10)
}

#### Discrimination - non-mechanistic####
for(task in names(task_colors)){
  beta_mean<-task_beta_prior[[task]][['mean']]
  beta_sd<-task_beta_prior[[task]][['sd']]

  discrimination_non_mechanistic<-
    tibble(idx=idx) %>% 
    rowwise() %>% 
    mutate(
      # Xi profile: mu1 is the level at condition 3 (AT=baseline), mu2/mu3 the independent outward
      # steps toward conditions 4/5, mu4/mu5 the outward steps toward conditions 2/1 (matches
      # discrimination_non_mechanistic.stan).
      mu1=rnorm(1,-2,0.5),
      mu2=rnorm(1,0,0.5),mu3=rnorm(1,0,0.5),mu4=rnorm(1,0,0.5),mu5=rnorm(1,0,0.5),
      # Beta profile: same chained construction, mu6 the condition-3 level
      mu6=rnorm(1,beta_mean,beta_sd),
      mu7=rnorm(1,0,0.5),mu8=rnorm(1,0,0.5),mu9=rnorm(1,0,0.5),mu10=rnorm(1,0,0.5),
      mu_logit_lambda=rnorm(1,-3.85+30*0.00171,0.72),
      mu_kappa=rnorm(1,0,0.5),
      # participant-level SDs, one per mu entry (tau ~ half-normal(0,1))
      tau1=abs(rnorm(1,0,0.5)),tau2=abs(rnorm(1,0,0.5)),tau3=abs(rnorm(1,0,0.5)),tau4=abs(rnorm(1,0,0.5)),tau5=abs(rnorm(1,0,0.5)),
      tau6=abs(rnorm(1,0,beta_sd)),tau7=abs(rnorm(1,0,beta_sd)),tau8=abs(rnorm(1,0,beta_sd)),tau9=abs(rnorm(1,0,beta_sd)),tau10=abs(rnorm(1,0,beta_sd)),
      tau_logit_lambda=abs(rnorm(1,0,0.72)),
      tau_kappa=abs(rnorm(1,0,0.5)),
      # group-level chained profiles (log scale of xi/beta)
      f_x3_g=mu1, f_x4_g=f_x3_g+mu2, f_x5_g=f_x4_g+mu3, f_x2_g=f_x3_g+mu4, f_x1_g=f_x2_g+mu5,
      f_b3_g=mu6, f_b4_g=f_b3_g+mu7, f_b5_g=f_b4_g+mu8, f_b2_g=f_b3_g+mu9, f_b1_g=f_b2_g+mu10,
      # one participant's chained draw (independent noise per step, ignoring cross-parameter correlation)
      f_x3_s=mu1+rnorm(1,0,tau1), f_x4_s=f_x3_s+mu2+rnorm(1,0,tau2), f_x5_s=f_x4_s+mu3+rnorm(1,0,tau3),
      f_x2_s=f_x3_s+mu4+rnorm(1,0,tau4), f_x1_s=f_x2_s+mu5+rnorm(1,0,tau5),
      f_b3_s=mu6+rnorm(1,0,tau6), f_b4_s=f_b3_s+mu7+rnorm(1,0,tau7), f_b5_s=f_b4_s+mu8+rnorm(1,0,tau8),
      f_b2_s=f_b3_s+mu9+rnorm(1,0,tau9), f_b1_s=f_b2_s+mu10+rnorm(1,0,tau10),
      lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
      kappa=rnorm(1,mu_kappa,tau_kappa)
    ) %>%
    ungroup() %>%
    full_join(grid_d) %>%
    mutate(
      c=at-29,
      log_xi_s=case_when(c==1~f_x1_s,c==2~f_x2_s,c==3~f_x3_s,c==4~f_x4_s,c==5~f_x5_s),
      log_beta_s =case_when(c==1~f_b1_s,c==2~f_b2_s,c==3~f_b3_s,c==4~f_b4_s,c==5~f_b5_s),
      log_xi_g=case_when(c==1~f_x1_g,c==2~f_x2_g,c==3~f_x3_g,c==4~f_x4_g,c==5~f_x5_g),
      log_beta_g =case_when(c==1~f_b1_g,c==2~f_b2_g,c==3~f_b3_g,c==4~f_b4_g,c==5~f_b5_g),

      interval_sign=sign(x),
      xi=exp(log_xi_s),
      beta=exp(log_beta_s),
      cx_s=abs(x)-xi,
      stim_rep=cx_s*inv_logit(cx_s*100),
      theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep+kappa),

      xi_g=exp(log_xi_g),
      beta_g=exp(log_beta_g),
      cx_g=abs(x)-xi_g,
      stim_rep_g=cx_g*inv_logit(cx_g*100),
      theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*beta_g*stim_rep_g+mu_kappa)
    ) %>%
    group_by(at,x) %>% 
    summarise(
      s_lb_ci_95=quantile(theta_s,0.025),
      s_ub_ci_95=quantile(theta_s,0.975),
      s_lb_ci_90=quantile(theta_s,0.05),
      s_ub_ci_90=quantile(theta_s,0.95),
      s_lb_ci_80=quantile(theta_s,0.1),
      s_ub_ci_80=quantile(theta_s,0.9),
      s_lb_ci_60=quantile(theta_s,0.20),
      s_ub_ci_60=quantile(theta_s,0.80),
      g_lb_ci_95=quantile(theta_g,0.025),
      g_ub_ci_95=quantile(theta_g,0.975),
      g_lb_ci_90=quantile(theta_g,0.05),
      g_ub_ci_90=quantile(theta_g,0.95),
      g_lb_ci_80=quantile(theta_g,0.1),
      g_ub_ci_80=quantile(theta_g,0.9),
      g_lb_ci_60=quantile(theta_g,0.20),
      g_ub_ci_60=quantile(theta_g,0.80),    
      s_med_ci_50 = quantile(theta_s, 0.5),
      g_med_ci_50 = quantile(theta_g, 0.5)
    ) %>% 
    pivot_longer(
      cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
      names_to = c("level","stat","ci"),
      names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
      values_to = "value"
    ) %>% 
    filter(TRUE) # response space: keep the full signed deviation range

  ribbons <- discrimination_non_mechanistic %>% filter(stat %in% c("lb","ub"))%>% 
    pivot_wider(
      names_from = "stat",
      values_from = "value"
    )
  medians <- discrimination_non_mechanistic %>% filter(stat == "med")

  p <- ribbons %>%
    ggplot()+
    geom_ribbon(
      aes(x=x,ymin = lb,ymax=ub,alpha=ci),
      fill=task_colors[[task]]
    )+
    geom_line(
      data=medians,
      aes(x=x,y=value),
      color=task_colors[[task]]
    )+
    geom_hline(
      aes(yintercept = .5),
      linetype='dotted',
      alpha=.5
    )+
    geom_hline(
      aes(yintercept = 1),
      linetype='dotted',
      alpha=.5
    )+
    geom_vline(aes(xintercept=0),linetype='dotted',alpha=.5)+
    scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
    scale_x_continuous(breaks=seq(-12,12,4))+
    scale_y_continuous(breaks=c(0,.5,1))+
    labs(
      alpha='',
      x='Signed target deviation (°C; + = deviating stimulus in 2nd interval)',
      y='P(chose 2nd interval)'
    )+
    theme_minimal()+
    facet_grid(rows=vars(level),cols=vars(at))

  ggsave(paste0('figures/d_n_ppp_',task,'.png'),plot=p,units = 'cm',width = 18,height = 10)
}

#### Rating - absolute coding####
rating_absolute<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_intercept=rnorm(1,-2,1),
    mu_log_slope=rnorm(1,-2,1),
    tau_intercept=abs(rnorm(1,0,1)),
    tau_log_slope=abs(rnorm(1,0,1)),
    intercept=rnorm(1,mu_intercept,tau_intercept),
    slope=exp(rnorm(1,mu_log_slope,tau_log_slope))
  ) %>%
  ungroup() %>%
  full_join(grid_r) %>%
  mutate(
    cx=(x-32),
    lr=intercept+slope*cx,
    theta_s=inv_logit(lr),
    lr_g=mu_intercept+exp(mu_log_slope)*cx,
    theta_g=inv_logit(lr_g)
  ) %>% 
  group_by(at,x) %>% 
  summarise(
    s_lb_ci_95=quantile(theta_s,0.025),
    s_ub_ci_95=quantile(theta_s,0.975),
    s_lb_ci_90=quantile(theta_s,0.05),
    s_ub_ci_90=quantile(theta_s,0.95),
    s_lb_ci_80=quantile(theta_s,0.1),
    s_ub_ci_80=quantile(theta_s,0.9),
    s_lb_ci_60=quantile(theta_s,0.20),
    s_ub_ci_60=quantile(theta_s,0.80),
    g_lb_ci_95=quantile(theta_g,0.025),
    g_ub_ci_95=quantile(theta_g,0.975),
    g_lb_ci_90=quantile(theta_g,0.05),
    g_ub_ci_90=quantile(theta_g,0.95),
    g_lb_ci_80=quantile(theta_g,0.1),
    g_ub_ci_80=quantile(theta_g,0.9),
    g_lb_ci_60=quantile(theta_g,0.20),
    g_ub_ci_60=quantile(theta_g,0.80),    
    s_med_ci_50 = quantile(theta_s, 0.5),
    g_med_ci_50 = quantile(theta_g, 0.5)
  ) %>% 
  pivot_longer(
    cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
    names_to = c("level","stat","ci"),
    names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
    values_to = "value"
  ) %>% 
  filter(x>at)

ribbons <- rating_absolute %>% filter(stat %in% c("lb","ub"))%>% 
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  )
medians <- rating_absolute %>% filter(stat == "med")

ribbons %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
  )+
  geom_line(
    data=medians,
    aes(x=x,y=value),
    color='#E69F00'
  )+
  geom_hline(
    aes(yintercept = 0),
    linetype='dotted',
    alpha=.5
  )+
  geom_hline(
    aes(yintercept = 1),
    linetype='dotted',
    alpha=.5
  )+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='Mean rating'
  )+
  theme_minimal()+
  facet_grid(rows=vars(level),cols=vars(at))

ggsave('figures/r_a_ppp.png',units = 'cm',width = 18,height = 10)

#### Rating - relative coding####
rating_relative<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_intercept=rnorm(1,-2,1),
    mu_log_slope=rnorm(1,-2,1),
    tau_intercept=abs(rnorm(1,0,1)),
    tau_log_slope=abs(rnorm(1,0,1)),
    intercept=rnorm(1,mu_intercept,tau_intercept),
    slope=exp(rnorm(1,mu_log_slope,tau_log_slope))
  ) %>%
  ungroup() %>%
  full_join(grid_r) %>%
  mutate(
    cx=(x-at),
    lr=intercept+slope*cx,
    theta_s=inv_logit(lr),
    
    lr_g=mu_intercept+exp(mu_log_slope)*cx,
    theta_g=inv_logit(lr_g)
  ) %>% 
  group_by(at,x) %>% 
  summarise(
    s_lb_ci_95=quantile(theta_s,0.025),
    s_ub_ci_95=quantile(theta_s,0.975),
    s_lb_ci_90=quantile(theta_s,0.05),
    s_ub_ci_90=quantile(theta_s,0.95),
    s_lb_ci_80=quantile(theta_s,0.1),
    s_ub_ci_80=quantile(theta_s,0.9),
    s_lb_ci_60=quantile(theta_s,0.20),
    s_ub_ci_60=quantile(theta_s,0.80),
    g_lb_ci_95=quantile(theta_g,0.025),
    g_ub_ci_95=quantile(theta_g,0.975),
    g_lb_ci_90=quantile(theta_g,0.05),
    g_ub_ci_90=quantile(theta_g,0.95),
    g_lb_ci_80=quantile(theta_g,0.1),
    g_ub_ci_80=quantile(theta_g,0.9),
    g_lb_ci_60=quantile(theta_g,0.20),
    g_ub_ci_60=quantile(theta_g,0.80),    
    s_med_ci_50 = quantile(theta_s, 0.5),
    g_med_ci_50 = quantile(theta_g, 0.5)
  ) %>% 
  pivot_longer(
    cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
    names_to = c("level","stat","ci"),
    names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
    values_to = "value"
  ) %>% 
  filter(x>at)

ribbons <- rating_relative %>% filter(stat %in% c("lb","ub"))%>% 
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  )
medians <- rating_relative %>% filter(stat == "med")

ribbons %>%
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
  )+
  geom_line(
    data=medians,
    aes(x=x,y=value),
    color='#E69F00'
  )+
  geom_hline(
    aes(yintercept = 0),
    linetype='dotted',
    alpha=.5
  )+
  geom_hline(
    aes(yintercept = 1),
    linetype='dotted',
    alpha=.5
  )+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='Mean rating'
  )+
  theme_minimal()+
  facet_grid(rows=vars(level),cols=vars(at))

ggsave('figures/r_r_ppp.png',units = 'cm',width = 18,height = 10)

#### Rating - non-mechanistic ####
rating_non_mechanistic<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    # Intercept profile: mu1 is the level at condition 3 (AT=baseline), mu2/mu3 the independent outward
    # steps toward conditions 4/5, mu4/mu5 the outward steps toward conditions 2/1 (matches
    # rating_non_mechanistic.stan).
    mu1=rnorm(1,-2,0.5),
    mu2=rnorm(1,0,0.5),mu3=rnorm(1,0,0.5),mu4=rnorm(1,0,0.5),mu5=rnorm(1,0,0.5),
    # Slope profile: same chained construction, mu6 the condition-3 (log-scale) level
    mu6=rnorm(1,-2,1),
    mu7=rnorm(1,0,0.5),mu8=rnorm(1,0,0.5),mu9=rnorm(1,0,0.5),mu10=rnorm(1,0,0.5),
    tau1=abs(rnorm(1,0,0.5)),tau2=abs(rnorm(1,0,0.5)),tau3=abs(rnorm(1,0,0.5)),tau4=abs(rnorm(1,0,0.5)),tau5=abs(rnorm(1,0,0.5)),
    tau6=abs(rnorm(1)),tau7=abs(rnorm(1)),tau8=abs(rnorm(1)),tau9=abs(rnorm(1)),tau10=abs(rnorm(1)),
    # group-level chained profiles
    f_i3_g=mu1, f_i4_g=f_i3_g+mu2, f_i5_g=f_i4_g+mu3, f_i2_g=f_i3_g+mu4, f_i1_g=f_i2_g+mu5,
    f_s3_g=mu6, f_s4_g=f_s3_g+mu7, f_s5_g=f_s4_g+mu8, f_s2_g=f_s3_g+mu9, f_s1_g=f_s2_g+mu10,
    # one participant's chained draw (independent noise per step, ignoring cross-parameter correlation)
    f_i3_s=mu1+rnorm(1,0,tau1), f_i4_s=f_i3_s+mu2+rnorm(1,0,tau2), f_i5_s=f_i4_s+mu3+rnorm(1,0,tau3),
    f_i2_s=f_i3_s+mu4+rnorm(1,0,tau4), f_i1_s=f_i2_s+mu5+rnorm(1,0,tau5),
    f_s3_s=mu6+rnorm(1,0,tau6), f_s4_s=f_s3_s+mu7+rnorm(1,0,tau7), f_s5_s=f_s4_s+mu8+rnorm(1,0,tau8),
    f_s2_s=f_s3_s+mu9+rnorm(1,0,tau9), f_s1_s=f_s2_s+mu10+rnorm(1,0,tau10)
  ) %>%
  ungroup() %>%
  full_join(grid_r) %>%
  mutate(
    c=at-29,
    intercept=case_when(c==1~f_i1_s,c==2~f_i2_s,c==3~f_i3_s,c==4~f_i4_s,c==5~f_i5_s),
    log_slope=case_when(c==1~f_s1_s,c==2~f_s2_s,c==3~f_s3_s,c==4~f_s4_s,c==5~f_s5_s),
    intercept_g=case_when(c==1~f_i1_g,c==2~f_i2_g,c==3~f_i3_g,c==4~f_i4_g,c==5~f_i5_g),
    log_slope_g=case_when(c==1~f_s1_g,c==2~f_s2_g,c==3~f_s3_g,c==4~f_s4_g,c==5~f_s5_g),

    slope=exp(log_slope),
    slope_g=exp(log_slope_g),

    cx=x-at,
    lr=intercept+slope*cx,
    theta_s=inv_logit(lr),

    lr_g=intercept_g+slope_g*cx,
    theta_g=inv_logit(lr_g)
  ) %>%
  group_by(at,x) %>% 
  summarise(
    s_lb_ci_95=quantile(theta_s,0.025),
    s_ub_ci_95=quantile(theta_s,0.975),
    s_lb_ci_90=quantile(theta_s,0.05),
    s_ub_ci_90=quantile(theta_s,0.95),
    s_lb_ci_80=quantile(theta_s,0.1),
    s_ub_ci_80=quantile(theta_s,0.9),
    s_lb_ci_60=quantile(theta_s,0.20),
    s_ub_ci_60=quantile(theta_s,0.80),
    g_lb_ci_95=quantile(theta_g,0.025),
    g_ub_ci_95=quantile(theta_g,0.975),
    g_lb_ci_90=quantile(theta_g,0.05),
    g_ub_ci_90=quantile(theta_g,0.95),
    g_lb_ci_80=quantile(theta_g,0.1),
    g_ub_ci_80=quantile(theta_g,0.9),
    g_lb_ci_60=quantile(theta_g,0.20),
    g_ub_ci_60=quantile(theta_g,0.80),    
    s_med_ci_50 = quantile(theta_s, 0.5),
    g_med_ci_50 = quantile(theta_g, 0.5)
  ) %>% 
  pivot_longer(
    cols = c(matches("(s|g)_(lb|ub|med)_ci_\\d+")),
    names_to = c("level","stat","ci"),
    names_pattern = "(s|g)_(lb|ub|med)_ci?_?(\\d+)?",
    values_to = "value"
  ) %>% 
  filter(x>at)

ribbons <- rating_non_mechanistic %>% filter(stat %in% c("lb","ub"))%>% 
  pivot_wider(
    names_from = "stat",
    values_from = "value"
  )
medians <- rating_non_mechanistic %>% filter(stat == "med")

ribbons%>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
  )+
  geom_line(
    data=medians,
    aes(x=x,y=value),
    color='#E69F00'
  )+
  geom_hline(
    aes(yintercept = 0),
    linetype='dotted',
    alpha=.5
  )+
  geom_hline(
    aes(yintercept = 1),
    linetype='dotted',
    alpha=.5
  )+
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(0,.5,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='Mean rating'
  )+
  theme_minimal()+
  facet_grid(rows=vars(level),cols=vars(at))

ggsave('figures/r_n_ppp.png',units = 'cm',width = 18,height = 10)
