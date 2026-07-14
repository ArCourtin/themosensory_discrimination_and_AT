# Script used to compare the three discrimination models and the three rating models
# (absolute, relative, non-mechanistic) via LOO, separately for the cold and warm adapting tasks.
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(loo)

rm(list=ls())

#### Functions ####
# results/loo/{domain}_{model}_{task}.rds: model 1=absolute, 2=relative, 3=non-mechanistic;
# task 1=cold, task 2=warm.
compare_models <- function(domain, task) {
  loos <- list(
    absolute = readRDS(paste0("results/loo/", domain, "_1_", task, ".rds")),
    relative = readRDS(paste0("results/loo/", domain, "_2_", task, ".rds")),
    `non-mechanistic` = readRDS(paste0("results/loo/", domain, "_3_", task, ".rds"))
  )
  loo_compare(loos) %>%
    as_tibble(rownames = "model") %>%
    mutate(z = elpd_diff / se_diff, p = pnorm(z))
}

#### Compare models per domain and task ####
domains <- c("discrimination", "rating")
tasks <- c(cold = 1, warm = 2)

comparisons <- map(domains, function(domain) map(tasks, ~compare_models(domain, .x)))
names(comparisons) <- domains

for (domain in domains) {
  for (task_name in names(comparisons[[domain]])) {
    cat("\n####", domain, "-", task_name, "task ####\n")
    print(comparisons[[domain]][[task_name]])
  }
}
