#This script is used to illustrate the priors used for the different models
#Author: A.S. Courtin
#Licence: MIT

library(tidyverse)

inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

M=10^3
x=seq(30,45,.1)
idx=1:M
at=30:34
grid<-expand.grid(x=x,idx=idx,at=at)

#### Discrimination - absolute coding####
discrimination_absolute<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_alpha=rnorm(1),
    mu_log_beta=rnorm(1),
    mu_logit_lambda=rnorm(1,-4,1),
    tau_alpha=abs(rnorm(1)),
    tau_log_beta=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    alpha=rnorm(1,mu_alpha,tau_alpha),
    beta=exp(rnorm(1,mu_log_beta,tau_log_beta)),
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2
  ) %>% 
  full_join(grid) %>% 
  mutate(
    cx_s=(x-34)-alpha,
    w_s=inv_logit(cx_s*100),
    theta_s=(1-w_s)*.5+(w_s-lambda)*pnorm(beta*cx_s),
    
    cx_g=(x-34)-mu_alpha,
    w_g=inv_logit(cx_g*100),
    theta_g=(1-w_g)*.5+(w_g-inv_logit(mu_logit_lambda)/2)*pnorm(exp(mu_log_beta)*cx_g)
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

discrimination_absolute %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(.5,.75,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='P(correct)'
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
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2
  ) %>% 
  full_join(grid) %>% 
  mutate(
    cx_s=(x-at)-alpha,
    w_s=inv_logit(cx_s*100),
    theta_s=(1-w_s)*.5+(w_s-lambda)*pnorm(beta*cx_s),
    
    cx_g=(x-at)-exp(mu_log_alpha),
    w_g=inv_logit(cx_g*100),
    theta_g=(1-w_g)*.5+(w_g-inv_logit(mu_logit_lambda)/2)*pnorm(exp(mu_log_beta)*cx_g)
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

discrimination_relative %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(.5,.75,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='P(correct)'
  )+
  theme_minimal()+
  facet_grid(rows=vars(level),cols=vars(at))

ggsave('figures/d_r_ppp.png',units = 'cm',width = 18,height = 10)

#### Discrimination - mixed coding####
discrimination_mixed<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_log_alpha=rnorm(1,-2,1),
    mu_log_sigma=rnorm(1,1,1),
    mu_logit_lambda=rnorm(1,-4,1),
    tau_log_alpha=abs(rnorm(1)),
    tau_log_sigma=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    alpha=exp(rnorm(1,mu_log_alpha,tau_log_alpha)),
    sigma=exp(rnorm(1,mu_log_sigma,tau_log_sigma)),
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2
  ) %>% 
  full_join(grid) %>% 
  mutate(
    beta=3/(2+sigma+alpha-(at-32)),
    cx_s=(x-at)-alpha,
    w_s=inv_logit(cx_s*100),
    theta_s=(1-w_s)*.5+(w_s-lambda)*pnorm(beta*cx_s),
    
    beta_g=3/(2+exp(mu_log_sigma)+exp(mu_log_alpha)-(at-32)),
    cx_g=(x-at)-exp(mu_log_alpha),
    w_g=inv_logit(cx_g*100),
    theta_g=(1-w_g)*.5+(w_g-inv_logit(mu_logit_lambda)/2)*pnorm(beta_g*cx_g)
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

discrimination_mixed %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(.5,.75,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='P(correct)'
  )+
  theme_minimal()+
  facet_grid(rows=vars(level),cols=vars(at))

ggsave('figures/d_m_ppp.png',units = 'cm',width = 18,height = 10)

#### Discrimination - non-mechanistic####
discrimination_non_mechanistic<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_alpha_int=rnorm(1),
    mu_alpha_d1=rnorm(1),
    mu_alpha_d2=rnorm(1),
    mu_alpha_d4=rnorm(1),
    mu_alpha_d5=rnorm(1),
    mu_log_beta_int=rnorm(1),
    mu_log_beta_d1=rnorm(1),
    mu_log_beta_d2=rnorm(1),
    mu_log_beta_d4=rnorm(1),
    mu_log_beta_d5=rnorm(1),
    mu_logit_lambda=rnorm(1,-4,1),
    tau_alpha_int=abs(rnorm(1)),
    tau_alpha_d1=abs(rnorm(1)),
    tau_alpha_d2=abs(rnorm(1)),
    tau_alpha_d4=abs(rnorm(1)),
    tau_alpha_d5=abs(rnorm(1)),
    tau_log_beta_int=abs(rnorm(1)),
    tau_log_beta_d1=abs(rnorm(1)),
    tau_log_beta_d2=abs(rnorm(1)),
    tau_log_beta_d4=abs(rnorm(1)),
    tau_log_beta_d5=abs(rnorm(1)),
    tau_logit_lambda=abs(rnorm(1)),
    
    alpha_int=rnorm(1,mu_alpha_int,tau_alpha_int),
    alpha_d1=rnorm(1,mu_alpha_d1,tau_alpha_d1),
    alpha_d2=rnorm(1,mu_alpha_d2,tau_alpha_d2),
    alpha_d4=rnorm(1,mu_alpha_d4,tau_alpha_d4),
    alpha_d5=rnorm(1,mu_alpha_d5,tau_alpha_d5),
    log_beta_int=rnorm(1,mu_log_beta_int,tau_log_beta_int),
    log_beta_d1=rnorm(1,mu_log_beta_d1,tau_log_beta_d1),
    log_beta_d2=rnorm(1,mu_log_beta_d2,tau_log_beta_d2),
    log_beta_d4=rnorm(1,mu_log_beta_d4,tau_log_beta_d4),
    log_beta_d5=rnorm(1,mu_log_beta_d5,tau_log_beta_d5),    
    lambda=inv_logit(rnorm(1,mu_logit_lambda,tau_logit_lambda))/2
  ) %>% 
  full_join(grid) %>% 
  mutate(
    at1=at==30,
    at2=at==31,
    at4=at==33,
    at5=at==34,
    
    alpha=alpha_int+at1*alpha_d1+at2*alpha_d2+at4*alpha_d4+at5*alpha_d5,
    beta=exp(log_beta_int+at1*log_beta_d1+at2*log_beta_d2+at4*log_beta_d4+at5*log_beta_d5),
    
    cx_s=(x-at)-alpha,
    w_s=inv_logit(cx_s*100),
    theta_s=(1-w_s)*.5+(w_s-lambda)*pnorm(beta*cx_s),
    
    alpha_g=mu_alpha_int+at1*mu_alpha_d1+at2*mu_alpha_d2+at4*mu_alpha_d4+at5*mu_alpha_d5,
    beta_g=exp(mu_log_beta_int+at1*mu_log_beta_d1+at2*mu_log_beta_d2+at4*mu_log_beta_d4+at5*mu_log_beta_d5),
    
    cx_g=(x-at)-alpha_g,
    w_g=inv_logit(cx_g*100),
    theta_g=(1-w_g)*.5+(w_g-inv_logit(mu_logit_lambda)/2)*pnorm(beta_g*cx_g)
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

discrimination_non_mechanistic %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
  scale_alpha_manual(labels=c('60% CI','80% CI','90% CI','95% CI'),values=c(.4,.3,.2,.1))+
  scale_x_continuous(breaks=seq(30,45,3))+
  scale_y_continuous(breaks=c(.5,.75,1))+
  labs(
    alpha='',
    x='Absolute target temperature',
    y='P(correct)'
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
    cx=(x-38),
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

rating_absolute %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
  full_join(grid) %>% 
  mutate(
    cx=(x-(at+8)),
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

rating_relative %>%
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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

#### Rating - mixed coding####
rating_mixed<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_log_alpha=rnorm(1,-2,1),
    mu_log_sigma=rnorm(1,1,1),
    tau_log_alpha=abs(rnorm(1)),
    tau_log_sigma=abs(rnorm(1)),
    alpha=exp(rnorm(1,mu_log_alpha,tau_log_alpha)),
    sigma=exp(rnorm(1,mu_log_sigma,tau_log_sigma)),
    mu_intercept=rnorm(1,-2,1),
    mu_log_slope=rnorm(1,-2,1),
    tau_intercept=abs(rnorm(1,0,1)),
    tau_log_slope=abs(rnorm(1,0,1)),
    intercept=rnorm(1,mu_intercept,tau_intercept),
    slope=exp(rnorm(1,mu_log_slope,tau_log_slope))
  ) %>% 
  full_join(grid) %>% 
  mutate(
    beta=3/(2+sigma+alpha-(at-32)),
    cx_s=x-(at+8)-alpha,
    beta_g=3/(2+exp(mu_log_sigma)+exp(mu_log_alpha)-(at-32)),
    cx_g=x-(at+8)-exp(mu_log_alpha),
    lr=intercept+slope*cx_s*beta,
    theta_s=inv_logit(lr),
    
    lr_g=mu_intercept+exp(mu_log_slope)*cx_g*beta_g,
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

rating_mixed %>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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

ggsave('figures/r_m_ppp.png',units = 'cm',width = 18,height = 10)

#### Rating - non-mechanistic ####
rating_non_mechanistic<-
  tibble(idx=idx) %>% 
  rowwise() %>% 
  mutate(
    mu_intercept_3=rnorm(1,-2,1),
    mu_intercept_d1=rnorm(1,0,1),
    mu_intercept_d2=rnorm(1,0,1),
    mu_intercept_d4=rnorm(1,0,1),
    mu_intercept_d5=rnorm(1,0,1),
    mu_log_slope_3=rnorm(1,-2,1),
    mu_log_slope_d1=rnorm(1,0,1),
    mu_log_slope_d2=rnorm(1,0,1),
    mu_log_slope_d4=rnorm(1,0,1),
    mu_log_slope_d5=rnorm(1,0,1),
    tau_intercept_3=abs(rnorm(1,0,1)),
    tau_intercept_d1=abs(rnorm(1,0,1)),
    tau_intercept_d2=abs(rnorm(1,0,1)),
    tau_intercept_d4=abs(rnorm(1,0,1)),
    tau_intercept_d5=abs(rnorm(1,0,1)),
    tau_log_slope_3=abs(rnorm(1,0,1)),
    tau_log_slope_d1=abs(rnorm(1,0,1)),
    tau_log_slope_d2=abs(rnorm(1,0,1)),
    tau_log_slope_d4=abs(rnorm(1,0,1)),
    tau_log_slope_d5=abs(rnorm(1,0,1)),
    intercept_3=rnorm(1,mu_intercept_3,tau_intercept_3),
    intercept_d1=rnorm(1,mu_intercept_d1,tau_intercept_d1),
    intercept_d2=rnorm(1,mu_intercept_d2,tau_intercept_d2),
    intercept_d4=rnorm(1,mu_intercept_d4,tau_intercept_d4),
    intercept_d5=rnorm(1,mu_intercept_d5,tau_intercept_d5),
    slope_3=rnorm(1,mu_log_slope_3,tau_log_slope_3),
    slope_d1=rnorm(1,mu_log_slope_d1,tau_log_slope_d1),
    slope_d2=rnorm(1,mu_log_slope_d2,tau_log_slope_d2),
    slope_d4=rnorm(1,mu_log_slope_d4,tau_log_slope_d4),
    slope_d5=rnorm(1,mu_log_slope_d5,tau_log_slope_d5)
  ) %>% 
  full_join(grid) %>% 
  mutate(
    at1=at==30,
    at2=at==31,
    at4=at==33,
    at5=at==34,
    
    intercept=intercept_3+at1*intercept_d1+at2*intercept_d2+at4*intercept_d4+at5*intercept_d5,
    slope=exp(slope_3+at1*slope_d1+at2*slope_d2+at4*slope_d4+at5*slope_d5),
    
    cx=x-(at+8),
    lr=intercept+slope*cx,
    theta_s=inv_logit(lr),
    
    intercept_g=mu_intercept_3+at1*mu_intercept_d1+at2*mu_intercept_d2+at4*mu_intercept_d4+at5*mu_intercept_d5,
    slope_g=exp(mu_log_slope_3+at1*mu_log_slope_d1+at2*mu_log_slope_d2+at4*mu_log_slope_d4+at5*mu_log_slope_d5),    
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
    g_ub_ci_60=quantile(theta_g,0.80)
  ) %>% 
  pivot_longer(
    cols = matches("(s|g)_(lb|ub)_ci_\\d+"),
    names_to = c("level","bound", "ci"), 
    names_pattern = "(s|g)_(lb|ub)_ci_(\\d+)", 
    values_to = "value"
  )  %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(level = recode(level,
                        "s" = "Participant",
                        "g" = "Group"))

rating_non_mechanistic%>% 
  filter(x>at) %>% 
  ggplot()+
  geom_ribbon(
    aes(x=x,ymin = lb,ymax=ub,alpha=ci),
    fill='#E69F00'
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
