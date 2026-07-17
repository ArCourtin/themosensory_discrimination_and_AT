# Script used to compare the four discrimination models and the four rating models (absolute
# with personal-baseline referencing, absolute with a fixed common reference, relative,
# non-mechanistic) via LOO, separately for the cold and warm adapting tasks. Models 1 and 2 share
# the same absolute stan model, differing only in whether recorded_baseline_temperature is each
# participant's measured baseline (1, fit by fit_{discrimination,rating}_models_SLURM.R) or a
# fixed value common to everyone (2, fit by fit_absolute_fixed_reference_models_SLURM.R).
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Edited with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(loo)

rm(list=ls())

#### Functions ####
# results/loo/{domain}_{model}_{task}.rds: model 1=absolute (personal baseline), 2=absolute
# (fixed reference), 3=relative, 4=non-mechanistic; task 1=cold, task 2=warm.
model_labels <- c(
  `1` = "absolute (personal baseline)",
  `2` = "absolute (fixed reference)",
  `3` = "relative",
  `4` = "non-mechanistic"
)

compare_models <- function(domain, task) {
  loos <- map(names(model_labels), function(m) {
    readRDS(paste0("results/loo/", domain, "_", m, "_", task, ".rds"))
  })
  names(loos) <- model_labels
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
