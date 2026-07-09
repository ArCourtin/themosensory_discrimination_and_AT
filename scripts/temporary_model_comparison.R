library(tidyverse)
library(cmdstanr)
library(loo)

l1 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/1_2.rds")
l2 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/2_2.rds")
l3 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/3_2.rds")

loo_compare(list(abs=l1,rel=l2,nm=l3))
lc_w<-
  loo_compare(list(abs=l1,rel=l2,nm=l3)) %>% 
  as.tibble() %>% 
  mutate(z=elpd_diff/se_diff,p=pnorm(z))
print(lc_w)


l1 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/1_1.rds")
l2 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/2_1.rds")
l3 <- readRDS("~/R/themosensory_discrimination_and_AT/results/loo/3_1.rds")

loo_compare(list(abs=l1,rel=l2,nm=l3))
lc_c<-
  loo_compare(list(abs=l1,rel=l2,nm=l3)) %>% 
  as.tibble() %>% 
  mutate(z=elpd_diff/se_diff,p=pnorm(z))
print(lc_c)

