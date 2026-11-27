#!/usr/bin/env Rscript
# ==============================================================================
# scripts/aggregate_results.R
# Merges Atomic .rds Checkpoints into Master Tables & Computes Jaccard/Nogueira Stability
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

source("R/05_symbolic_canonicalize.R")
source("R/06_metrics_consistency.R")
source("R/01_ode_simulators.R")
source("R/10_real_data_experiments.R")

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

  # Enforce audited canonicalization + failed-run NA scoring across all rows
  reg <- if (suite == "synth") get_systems_registry() else get_real_datasets_registry()
  for (i in seq_len(nrow(dt))) {
    st <- if ("Status" %in% names(dt)) dt$Status[i] else "OK"
    if (is.na(st) || st != "OK") {
      dt[i, `:=`(
        TPR = 0.0,
        FDR = NA_real_,
        F1 = 0.0,
        Exact_Recovery = 0L,
        SHD = NA_integer_,
        OC = NA_integer_,
        n_active_terms = NA_integer_
      )]
      next
    }
    key <- if (suite == "synth") as.character(dt$System[i]) else as.character(dt$Dataset[i])
    if (!key %in% names(reg)) next
    gt_list <- reg[[key]]$ground_truth
    m_eq <- length(gt_list)
    raw_str <- dt$Discovered_Atoms[i]
    disc_list <- replicate(m_eq, character(0), simplify = FALSE)
    if (!is.na(raw_str) && nzchar(raw_str)) {
      toks <- strsplit(raw_str, ";", fixed = TRUE)[[1]]
      for (tk in toks) {
        if (grepl("^eq([0-9]+)::(.*)$", tk)) {
          eq_i <- as.integer(sub("^eq([0-9]+)::(.*)$", "\\1", tk))
          atom_i <- sub("^eq([0-9]+)::(.*)$", "\\2", tk)
          if (eq_i >= 1L && eq_i <= m_eq) {
            disc_list[[eq_i]] <- c(disc_list[[eq_i]], atom_i)
          }
        }
      }
    }
    m_res <- compute_structure_metrics(disc_list, gt_list)
    dt[i, `:=`(
      TPR = m_res$TPR,
      FDR = m_res$FDR,
      F1 = m_res$F1,
      Exact_Recovery = m_res$Exact_Recovery,
      SHD = m_res$SHD,
      OC = m_res$OC,
      n_active_terms = m_res$n_active_terms,
      Discovered_Atoms = paste(m_res$S_hat, collapse = ";")
    )]
  }

  out_rds <- file.path("results/aggregated", paste0(suite, "_results.rds"))
  out_csv <- file.path("results/aggregated", paste0(suite, "_results.csv"))
  saveRDS(dt, out_rds)
  fwrite(dt, out_csv)
  message("[OK] Aggregated ", nrow(dt), " tasks for suite '", suite, "' -> ", out_csv)

  # Compute selection stability (Jaccard & Nogueira) across successful replications per condition x Method
  split_sets <- function(s_vec) {
    lapply(s_vec, function(s) {
      if (is.na(s) || !nzchar(s)) character(0) else strsplit(s, ";", fixed = TRUE)[[1]]
    })
  }

  if (suite == "synth") {
    stab_dt <- dt[, {
      ok_idx <- which(Status == "OK")
      sets <- split_sets(Discovered_Atoms[ok_idx])
      nog <- compute_nogueira_stability(sets, d_total = 100L)
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
      ok_idx <- which(Status == "OK")
      sets <- split_sets(Discovered_Atoms[ok_idx])
      nog <- compute_nogueira_stability(sets, d_total = 100L)
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

args_cli <- commandArgs(trailingOnly = TRUE)
s_idx <- which(args_cli == "--suite")
suite_req <- if (length(s_idx) > 0L && length(args_cli) > s_idx[1]) args_cli[s_idx[1] + 1L] else "both"

if (suite_req %in% c("synth", "both")) aggregate_suite("synth")
if (suite_req %in% c("real", "both")) aggregate_suite("real")

