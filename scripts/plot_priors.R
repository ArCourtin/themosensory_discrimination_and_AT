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

#### Discrimination - absolute coding####
discrimination_absolute<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_log_alpha=rnorm(1,-2,1),
    mu_rho=rnorm(1,0,2),
    mu_log_beta=rnorm(1),
    mu_logit_lambda=rnorm(1,-4,1),
    tau_log_alpha=abs(rnorm(1,0,1)),
    tau_rho=abs(rnorm(1,0,2)),
    tau_log_beta=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    alpha=exp(rnorm(1,mu_log_alpha,tau_log_alpha)),
    rho=rnorm(1,mu_rho,tau_rho),
    beta=exp(rnorm(1,mu_log_beta,tau_log_beta)),
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
    mu_kappa=rnorm(1,0,1),
    tau_kappa=abs(rnorm(1,0,1)),
    kappa=rnorm(1,mu_kappa,tau_kappa)
  ) %>%
  full_join(grid_d) %>%
  mutate(
    interval_sign=sign(x),
    target=at+abs(x),
    cx_s=target-32-rho,
    stim_rep=cx_s*inv_logit(cx_s*100),
    adapt_stim_rep=stim_rep*inv_logit(100*(stim_rep-(at-32+alpha-rho))),
    theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*adapt_stim_rep-kappa),
    
    cx_g=target-32-mu_rho,
    stim_rep_g=cx_g*inv_logit(cx_g*100),
    adapt_stim_rep_g=stim_rep_g*inv_logit(100*(stim_rep_g-(at-32+exp(mu_log_alpha)-mu_rho))),
    theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*exp(mu_log_beta)*adapt_stim_rep_g-mu_kappa)
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

ggsave('figures/d_a_ppp.png',units = 'cm',width = 18,height = 10)

#### Discrimination - relative coding####
discrimination_relative<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_log_alpha=rnorm(1,-2,1),
    mu_log_beta=rnorm(1),
    mu_logit_lambda=rnorm(1,-4,1),
    tau_log_alpha=abs(rnorm(1)),
    tau_log_beta=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    alpha=exp(rnorm(1,mu_log_alpha,tau_log_alpha)),
    beta=exp(rnorm(1,mu_log_beta,tau_log_beta)),
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
    mu_kappa=rnorm(1,0,1),
    tau_kappa=abs(rnorm(1,0,1)),
    kappa=rnorm(1,mu_kappa,tau_kappa)
  ) %>%
  full_join(grid_d) %>%
  mutate(
    interval_sign=sign(x),
    cx_s=abs(x)-alpha,
    stim_rep=cx_s*inv_logit(cx_s*100),
    theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa),

    cx_g=abs(x)-exp(mu_log_alpha),
    stim_rep_g=cx_g*inv_logit(cx_g*100),
    theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*exp(mu_log_beta)*stim_rep_g-mu_kappa)
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

ggsave('figures/d_r_ppp.png',units = 'cm',width = 18,height = 10)

#### Discrimination - non-mechanistic####
discrimination_non_mechanistic<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    # RW2 group profile of log-alpha across the ordered AT axis (2 anchors + curvature penalised by sigma_a)
    a0=rnorm(1,-2,1),
    a1=rnorm(1,0,0.5),
    a_c1=rnorm(1),a_c2=rnorm(1),a_c3=rnorm(1),
    sigma_a=abs(rnorm(1,0,0.3)),
    # RW2 group profile of log-beta
    b0=rnorm(1,0,1),
    b1=rnorm(1,0,0.5),
    b_c1=rnorm(1),b_c2=rnorm(1),b_c3=rnorm(1),
    sigma_b=abs(rnorm(1,0,0.3)),
    mu_logit_lambda=rnorm(1,-4,1),
    mu_kappa=rnorm(1,0,1),
    # per-condition participant SDs, plus lambda/kappa SDs (tau ~ half-normal(0,1))
    tau_a1=abs(rnorm(1)),tau_a2=abs(rnorm(1)),tau_a3=abs(rnorm(1)),tau_a4=abs(rnorm(1)),tau_a5=abs(rnorm(1)),
    tau_b1=abs(rnorm(1)),tau_b2=abs(rnorm(1)),tau_b3=abs(rnorm(1)),tau_b4=abs(rnorm(1)),tau_b5=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    tau_kappa=abs(rnorm(1)),
    # group-level RW2 profiles (log scale of alpha/beta)
    f_a1=a0, f_a2=a0+a1,
    f_a3=2*f_a2-f_a1+sigma_a*a_c1,
    f_a4=2*f_a3-f_a2+sigma_a*a_c2,
    f_a5=2*f_a4-f_a3+sigma_a*a_c3,
    f_b1=b0, f_b2=b0+b1,
    f_b3=2*f_b2-f_b1+sigma_b*b_c1,
    f_b4=2*f_b3-f_b2+sigma_b*b_c2,
    f_b5=2*f_b4-f_b3+sigma_b*b_c3,
    # one participant draw per condition (log scale)
    la1=rnorm(1,f_a1,tau_a1),la2=rnorm(1,f_a2,tau_a2),la3=rnorm(1,f_a3,tau_a3),la4=rnorm(1,f_a4,tau_a4),la5=rnorm(1,f_a5,tau_a5),
    lb1=rnorm(1,f_b1,tau_b1),lb2=rnorm(1,f_b2,tau_b2),lb3=rnorm(1,f_b3,tau_b3),lb4=rnorm(1,f_b4,tau_b4),lb5=rnorm(1,f_b5,tau_b5),
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2,
    kappa=rnorm(1,mu_kappa,tau_kappa)
  ) %>%
  full_join(grid_d) %>%
  mutate(
    c=at-29,
    log_alpha_s=case_when(c==1~la1,c==2~la2,c==3~la3,c==4~la4,c==5~la5),
    log_beta_s =case_when(c==1~lb1,c==2~lb2,c==3~lb3,c==4~lb4,c==5~lb5),
    log_alpha_g=case_when(c==1~f_a1,c==2~f_a2,c==3~f_a3,c==4~f_a4,c==5~f_a5),
    log_beta_g =case_when(c==1~f_b1,c==2~f_b2,c==3~f_b3,c==4~f_b4,c==5~f_b5),

    interval_sign=sign(x),
    alpha=exp(log_alpha_s),
    beta=exp(log_beta_s),
    cx_s=abs(x)-alpha,
    stim_rep=cx_s*inv_logit(cx_s*100),
    theta_s=lambda+(1-2*lambda)*pnorm(interval_sign*beta*stim_rep-kappa),

    alpha_g=exp(log_alpha_g),
    beta_g=exp(log_beta_g),
    cx_g=abs(x)-alpha_g,
    stim_rep_g=cx_g*inv_logit(cx_g*100),
    theta_g=inv_logit(mu_logit_lambda)/2+(1-2*inv_logit(mu_logit_lambda)/2)*pnorm(interval_sign*beta_g*stim_rep_g-mu_kappa)
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

ggsave('figures/d_n_ppp.png',units = 'cm',width = 18,height = 10)

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
  full_join(grid) %>% 
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
    # RW2 group profile of the latent intercept across the ordered AT axis (identity scale)
    i0=rnorm(1,-2,1),
    i1=rnorm(1,0,0.5),
    i_c1=rnorm(1),i_c2=rnorm(1),i_c3=rnorm(1),
    sigma_i=abs(rnorm(1,0,0.3)),
    # RW2 group profile of the latent log-slope
    s0=rnorm(1,-2,1),
    s1=rnorm(1,0,0.5),
    s_c1=rnorm(1),s_c2=rnorm(1),s_c3=rnorm(1),
    sigma_s=abs(rnorm(1,0,0.3)),
    tau_i1=abs(rnorm(1)),tau_i2=abs(rnorm(1)),tau_i3=abs(rnorm(1)),tau_i4=abs(rnorm(1)),tau_i5=abs(rnorm(1)),
    tau_s1=abs(rnorm(1)),tau_s2=abs(rnorm(1)),tau_s3=abs(rnorm(1)),tau_s4=abs(rnorm(1)),tau_s5=abs(rnorm(1)),
    # group-level RW2 profiles
    f_i1=i0, f_i2=i0+i1,
    f_i3=2*f_i2-f_i1+sigma_i*i_c1,
    f_i4=2*f_i3-f_i2+sigma_i*i_c2,
    f_i5=2*f_i4-f_i3+sigma_i*i_c3,
    f_s1=s0, f_s2=s0+s1,
    f_s3=2*f_s2-f_s1+sigma_s*s_c1,
    f_s4=2*f_s3-f_s2+sigma_s*s_c2,
    f_s5=2*f_s4-f_s3+sigma_s*s_c3,
    # one participant draw per condition
    ii1=rnorm(1,f_i1,tau_i1),ii2=rnorm(1,f_i2,tau_i2),ii3=rnorm(1,f_i3,tau_i3),ii4=rnorm(1,f_i4,tau_i4),ii5=rnorm(1,f_i5,tau_i5),
    ss1=rnorm(1,f_s1,tau_s1),ss2=rnorm(1,f_s2,tau_s2),ss3=rnorm(1,f_s3,tau_s3),ss4=rnorm(1,f_s4,tau_s4),ss5=rnorm(1,f_s5,tau_s5)
  ) %>%
  full_join(grid_r) %>% 
  mutate(
    c=at-29,
    intercept=case_when(c==1~ii1,c==2~ii2,c==3~ii3,c==4~ii4,c==5~ii5),
    log_slope=case_when(c==1~ss1,c==2~ss2,c==3~ss3,c==4~ss4,c==5~ss5),
    intercept_g=case_when(c==1~f_i1,c==2~f_i2,c==3~f_i3,c==4~f_i4,c==5~f_i5),
    log_slope_g=case_when(c==1~f_s1,c==2~f_s2,c==3~f_s3,c==4~f_s4,c==5~f_s5),

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
