#!/usr/bin/env Rscript
# ==============================================================================
# scripts/aggregate_results.R
# Merges Atomic .rds Checkpoints into Master Tables & Computes Jaccard/Nogueira Stability
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

source("R/06_metrics_consistency.R")

dir.create("results/aggregated", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

aggregate_suite <- function(suite = c("synth", "real")) {
  suite <- match.arg(suite)
  ckpt_dir <- paste0("results/checkpoints_", suite)
  files <- list.files(ckpt_dir, pattern = "^task_.*\\.rds$", full.names = TRUE)
  if (length(files) == 0L) {
    message("[AGGREGATE] No checkpoints found in ", ckpt_dir)
    return(NULL)
  }
  rows_list <- lapply(files, function(f) {
    tryCatch(readRDS(f), error = function(e) {
      warning("[AGGREGATE] Skipping corrupted checkpoint: ", f, " (", e$message, ")")
      NULL
    })
  })
  rows_list <- rows_list[!vapply(rows_list, is.null, logical(1))]
  if (length(rows_list) == 0L) {
    message("[AGGREGATE] No valid checkpoints found in ", ckpt_dir)
    return(NULL)
  }
  dt <- rbindlist(rows_list, fill = TRUE)
  setorder(dt, task_id)

  out_rds <- file.path("results/aggregated", paste0(suite, "_results.rds"))
  out_csv <- file.path("results/aggregated", paste0(suite, "_results.csv"))
  saveRDS(dt, out_rds)
  fwrite(dt, out_csv)
  message("[OK] Aggregated ", nrow(dt), " tasks for suite '", suite, "' -> ", out_csv)

  # Compute selection stability (Jaccard & Nogueira) across seeds/folds per condition x Method
  split_sets <- function(s_vec) {
    lapply(s_vec, function(s) {
      if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, ";", fixed = TRUE)[[1]]
    })
  }

  if (suite == "synth") {
    stab_dt <- dt[, {
      sets <- split_sets(Discovered_Atoms)
      nog <- compute_nogueira_stability(sets)
      .(
        n_reps = .N,
        Mean_TPR = mean(TPR, na.rm = TRUE),
        Mean_FDR = mean(FDR, na.rm = TRUE),
        Mean_F1 = mean(F1, na.rm = TRUE),
        Exact_Rate = mean(Exact_Recovery, na.rm = TRUE),
        Mean_SHD = mean(SHD, na.rm = TRUE),
        Mean_OC = mean(OC, na.rm = TRUE),
        Mean_PIP_Brier = mean(PIP_Brier, na.rm = TRUE),
        Mean_R2_in = mean(R2_test_in, na.rm = TRUE),
        Mean_R2_ood = mean(R2_test_ood, na.rm = TRUE),
        Mean_R2_traj = mean(R2_traj, na.rm = TRUE),
        Mean_T_valid = mean(T_valid, na.rm = TRUE),
        Mean_Coverage_95 = mean(Coverage_95_ood, na.rm = TRUE),
        Jaccard_Stability = compute_jaccard_stability(sets),
        Nogueira_Phi = nog$phi,
        Nogueira_CI_Lower = nog$ci_lower,
        Nogueira_CI_Upper = nog$ci_upper,
        Mean_Runtime_s = mean(Runtime_s, na.rm = TRUE)
      )
    }, by = .(System, SystemGroup, SampleSize, Noise, Method)]
    fwrite(stab_dt, "results/aggregated/synth_cell_summary_stability.csv")
  } else {
    stab_dt <- dt[, {
      sets <- split_sets(Discovered_Atoms)
      nog <- compute_nogueira_stability(sets)
      .(
        n_folds = .N,
        Mean_TPR = mean(TPR, na.rm = TRUE),
        Mean_FDR = mean(FDR, na.rm = TRUE),
        Mean_F1 = mean(F1, na.rm = TRUE),
        Mean_OC = mean(OC, na.rm = TRUE),
        Mean_R2_ood = mean(R2_test_ood, na.rm = TRUE),
        Mean_R2_traj = mean(R2_traj, na.rm = TRUE),
        Mean_NRMSE_traj = mean(NRMSE_traj, na.rm = TRUE),
        Mean_Coverage_95 = mean(Coverage_95_ood, na.rm = TRUE),
        Jaccard_Stability = compute_jaccard_stability(sets),
        Nogueira_Phi = nog$phi,
        Nogueira_CI_Lower = nog$ci_lower,
        Nogueira_CI_Upper = nog$ci_upper,
        Mean_Runtime_s = mean(Runtime_s, na.rm = TRUE)
      )
    }, by = .(Dataset, Domain, DiffMethod, Method)]
    fwrite(stab_dt, "results/aggregated/real_cell_summary_stability.csv")
  }
  dt
}

aggregate_suite("synth")
aggregate_suite("real")

