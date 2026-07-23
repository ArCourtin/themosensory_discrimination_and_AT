# Per-participant decomposition of the discrimination/rating LOO model comparisons run in
# model_comparison.R. loo_compare's elpd_diff pools every trial from every participant together,
# so a single participant (e.g. one who engaged poorly with the task) contributing many trials,
# or trials on which the models disagree sharply, can dominate the pooled comparison even if most
# participants individually favour the other model. This script instead: (1) splits each model's
# pointwise elpd_loo by participant to show how many participants individually favour each model,
# and (2) reports a participant-clustered SE for elpd_diff (sqrt(P_participants * var(per-
# participant diff))) alongside the usual trial-clustered SE (sqrt(N_trials * var(per-trial
# diff))) that loo_compare itself reports. The participant-clustered SE doesn't assume
# independence between trials sharing the same participant's hierarchical parameters, which the
# trial-clustered SE does.
#
# Trial counts per participant vary a lot post-exclusion, but not as random per-trial noise: the
# adaptation response was only saved for round 1, so an entire AT-condition block (~30 trials) is
# dropped whenever adaptation wasn't complete straight away, participant by participant. A
# participant left with only 1 surviving block has data from a single AT condition, and the three
# competing models differ specifically in how they predict behaviour *changes across* AT
# conditions — with no within-participant AT variation left, such a participant is close to
# uninformative for distinguishing the models, not just noisier. min_blocks below excludes
# participants with fewer surviving blocks than that from the win-count/SE, reporting the
# excluded count for transparency; the full per-participant table (with a block-count column) is
# still printed for all participants regardless.
# Author: Arthur S. Courtin
# License: MIT (see LICENSE file)
# Written with the assistance of Claude Code (Anthropic).

#### Set-up environment ####
library(tidyverse)
library(loo)

rm(list=ls())

min_blocks <- 2

#### Functions ####
# results/loo/{domain}_{model}_{task}.rds: model 1=absolute (personal baseline), 2=absolute
# (fixed reference), 3=relative, 4=non-mechanistic; task 1=cold, task 2=warm.
model_labels <- c(
  `1` = "absolute (personal baseline)",
  `2` = "absolute (fixed reference)",
  `3` = "relative",
  `4` = "non-mechanistic"
)

# Reproduces the per-trial participant index and AT-condition block feeding the `participant` /
# `adapting_temperature_idx` entries of the Stan data list in fit_discrimination_models_SLURM.R /
# fit_rating_models_SLURM.R / fit_absolute_fixed_reference_models_SLURM.R (same filters, same
# renumbering, same row order), so each row of a model's pointwise elpd_loo can be attributed
# back to a participant and an AT-condition block. `task_code` here is the saved-filename task
# suffix (1 = cold, 2 = warm); the raw CSV's own `task` column is coded 0 = cold, 1 = warm
# (task_code - 1), matching the fitting scripts' `t`.
load_trial_metadata <- function(domain, task_code) {
  if (domain == "discrimination") {
    data <- read_csv("data/d_at_2ifc_af.csv", show_col_types = FALSE) %>%
      filter(baseline_flag == 0, deviation_flag == 0)
  } else {
    data <- read_csv("data/d_at_ratings_af.csv", show_col_types = FALSE) %>%
      filter(baseline_flag == 0, deviation_flag == 0, confirmed == 1, !is.na(recorded_temperature))
  }
  participant_ids <- unique(data$participant)
  for (pdx in seq_along(participant_ids)) {
    data$participant[data$participant == participant_ids[pdx]] <- pdx
  }
  data %>%
    filter(task == task_code - 1) %>%
    mutate(adapting_temperature_idx = 3 + round(adapting - baseline)) %>%
    select(participant, adapting_temperature_idx)
}

load_pointwise_elpd <- function(domain, model, task, n_expected) {
  loo_obj <- readRDS(paste0("results/loo/", domain, "_", model, "_", task, ".rds"))
  elpd <- loo_obj$pointwise[, "elpd_loo"]
  stopifnot(length(elpd) == n_expected)
  elpd
}

# elpd_diff and se_diff_trial_clustered reproduce what loo_compare itself would report for this
# pair (trials treated as the independent unit). se_diff_participant_clustered instead treats each
# participant's summed diff as the independent unit; z/p are computed from the latter.
compare_pair <- function(elpd_a, elpd_b, participant_idx, label_a, label_b) {
  diff_trial <- elpd_a - elpd_b
  n_trial <- length(diff_trial)
  diff_participant <- tapply(diff_trial, participant_idx, sum)
  p_participants <- length(diff_participant)

  elpd_diff <- sum(diff_trial)
  se_participant_clustered <- sqrt(p_participants * var(diff_participant))

  tibble(
    model_a = label_a,
    model_b = label_b,
    elpd_diff = elpd_diff,
    se_diff_trial_clustered = sqrt(n_trial * var(diff_trial)),
    se_diff_participant_clustered = se_participant_clustered,
    n_favouring_a = sum(diff_participant > 0),
    n_favouring_b = sum(diff_participant < 0),
    z_participant_clustered = elpd_diff / se_participant_clustered,
    p_value = pnorm(elpd_diff / se_participant_clustered)
  )
}

#### Compare models per domain and task ####
# Run twice: once across all four models, once with the non-mechanistic model (4) left out, same
# as model_comparison.R, since it's a flexible envelope rather than a competing mechanistic
# account and can dominate comparisons against the three mechanistic models by fit alone.
domains <- c("discrimination", "rating")
tasks <- c(cold = 1, warm = 2)
model_sets <- list(
  `with non-mechanistic` = names(model_labels),
  `without non-mechanistic` = setdiff(names(model_labels), "4")
)

for (set_name in names(model_sets)) {
  models <- model_sets[[set_name]]
  for (domain in domains) {
    for (task_name in names(tasks)) {
      task <- tasks[[task_name]]
      trial_meta <- load_trial_metadata(domain, task)
      participant_idx <- trial_meta$participant
      n_trial <- nrow(trial_meta)

      elpd_list <- map(models, ~ load_pointwise_elpd(domain, .x, task, n_trial))
      names(elpd_list) <- model_labels[models]

      participant_ids_sorted <- sort(unique(participant_idx))
      n_blocks <- as.numeric(tapply(trial_meta$adapting_temperature_idx, participant_idx, n_distinct))
      n_trials_participant <- as.numeric(table(participant_idx))
      participant_sums <- map_dfc(elpd_list, ~ as.numeric(tapply(.x, participant_idx, sum)))
      participant_sums <- bind_cols(
        participant = participant_ids_sorted,
        n_blocks = n_blocks,
        n_trials = n_trials_participant,
        participant_sums
      )

      cat("\n####", set_name, "-", domain, "-", task_name, "task ####\n")
      cat("Per-participant summed elpd_loo by model (n_blocks = surviving AT-condition blocks, out of 5):\n")
      print(participant_sums)

      elpd_cols <- setdiff(names(participant_sums), c("participant", "n_blocks", "n_trials"))
      eligible <- participant_sums %>% filter(n_blocks >= min_blocks)
      n_excluded <- nrow(participant_sums) - nrow(eligible)
      if (n_excluded > 0) {
        cat("\n", n_excluded, "participant(s) with fewer than", min_blocks,
            "surviving blocks excluded from the win-count/SE below (shown above for reference).\n")
      }

      best_model <- apply(
        as.matrix(select(eligible, all_of(elpd_cols))), 1,
        function(row) names(row)[which.max(row)]
      )
      cat("\nParticipants (n_blocks >=", min_blocks, ") individually best explained by each model:\n")
      print(table(best_model))

      pooled_ranking <- sort(colSums(select(eligible, all_of(elpd_cols))), decreasing = TRUE)
      top_two <- names(pooled_ranking)[1:2]
      keep_trial <- participant_idx %in% eligible$participant

      cat("\nClustered comparison (n_blocks >=", min_blocks, "), top two by pooled elpd_loo (",
          top_two[1], "vs", top_two[2], "):\n")
      print(compare_pair(
        elpd_list[[top_two[1]]][keep_trial],
        elpd_list[[top_two[2]]][keep_trial],
        participant_idx[keep_trial],
        top_two[1], top_two[2]
      ))
    }
  }
}
