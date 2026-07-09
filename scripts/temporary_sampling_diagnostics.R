library(cmdstanr)
library(bayesplot)

fit <- readRDS("~/R/themosensory_discrimination_and_AT/results/fits/3_2.rds")
mcmc_intervals(fit$draws('alpha'))
mcmc_intervals(fit$draws('beta'))
mcmc_intervals(fit$draws('lambda'))



np <- nuts_params(fit)
mcmc_pairs(fit$draws(c('mu[1]','tau[1]')),np=np,transformations = list('tau[1]'='log'))
mcmc_pairs(fit$draws(c('mu[2]','tau[2]')),np=np,transformations = list('tau[2]'='log'))
mcmc_pairs(fit$draws(c('mu[3]','tau[3]')),np=np,transformations = list('tau[3]'='log'))
mcmc_pairs(fit$draws(c('mu[4]','tau[4]')),np=np,transformations = list('tau[4]'='log'))
mcmc_pairs(fit$draws(c('mu[5]','tau[5]')),np=np)
mcmc_pairs(fit$draws(c('mu[1]','tau[1]')),np=np)
