# Script used to compare the three discrimination models and the three rating models (absolute,
# relative, non-mechanistic) via LOO, separately for the cold and warm adapting tasks. Each
# domain/task comparison is run twice: once across all three models, once with the non-mechanistic
# model (3) excluded, since it's a flexible envelope rather than a competing mechanistic account.
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
model_labels <- c(
  `1` = "absolute",
  `2` = "relative",
  `3` = "non-mechanistic"
)

compare_models <- function(domain, task, models = names(model_labels)) {
  loos <- map(models, function(m) {
    readRDS(paste0("results/loo/", domain, "_", m, "_", task, ".rds"))
  })
  names(loos) <- model_labels[models]
  loo_compare(loos) %>%
    as_tibble(rownames = "model") %>%
    mutate(z = elpd_diff / se_diff, p = pnorm(z))
}

#### Compare models per domain and task ####
# Run twice: once across all three models, once with the non-mechanistic model (3) left out,
# since it's a flexible envelope rather than a competing mechanistic account and can dominate
# LOO comparisons against the two mechanistic models by fit alone.
domains <- c("discrimination", "rating")
tasks <- c(cold = 1, warm = 2)
model_sets <- list(
  `with non-mechanistic` = names(model_labels),
  `without non-mechanistic` = setdiff(names(model_labels), "3")
)

for (set_name in names(model_sets)) {
  models <- model_sets[[set_name]]
  for (domain in domains) {
    for (task_name in names(tasks)) {
      cat("\n####", set_name, "-", domain, "-", task_name, "task ####\n")
      print(compare_models(domain, tasks[[task_name]], models))
    }
  }
}
