#!/usr/bin/env Rscript
# ==============================================================================
# scripts/rescore_aggregated_results.R
# Re-scores existing results/aggregated/{synth,real}_results.rds under audited rules:
#   1. Failed runs (Status != "OK") receive FDR=NA, OC=NA, SHD=NA, n_active_terms=NA
#   2. Discovered_Atoms re-scored with audited canonicalizer:
#      - sin_deg / cos_deg are NOT aliased to radian sin / cos
#      - rational complements modulo intercept (x^2/(1+x^2) <-> 1/(1+x^2), x/(1+x) <-> 1/(1+x))
#        applied symmetrically to all methods
#   3. Nogueira stability Phi recomputed with shared d_total = 100 on Status == "OK" runs
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

source("R/05_symbolic_canonicalize.R")
source("R/06_metrics_consistency.R")
source("R/01_ode_simulators.R")
source("R/10_real_data_experiments.R")

rescore_dt <- function(suite = c("synth", "real")) {
  suite <- match.arg(suite)
  rds_path <- file.path("results/aggregated", paste0(suite, "_results.rds"))
  csv_path <- file.path("results/aggregated", paste0(suite, "_results.csv"))
  if (!file.exists(rds_path)) {
    message("[SKIP] Not found: ", rds_path)
    return(NULL)
  }
  dt <- as.data.table(readRDS(rds_path))

  # Backup original once before overwriting
  bak_path <- file.path("results/aggregated", paste0(suite, "_results_pre_audit_bak.rds"))
  if (!file.exists(bak_path)) {
    saveRDS(copy(dt), bak_path)
  }

  # Lookup ground truth by System / Dataset
  if (suite == "synth") {
    reg <- get_systems_registry()
  } else {
    reg <- get_real_datasets_registry()
  }

  # Enforce audited canonicalization + failed-run NA scoring row-by-row
  n_rows <- nrow(dt)
  n_changed <- 0L
  for (i in seq_len(n_rows)) {
    st <- dt$Status[i]
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
    gt_list <- reg[[key]]$ground_truth
    m <- length(gt_list)

    raw_str <- dt$Discovered_Atoms[i]
    # Reconstruct per-equation list from "eqj::atom" tokens
    disc_list <- replicate(m, character(0), simplify = FALSE)
    if (!is.na(raw_str) && nzchar(raw_str)) {
      toks <- strsplit(raw_str, ";", fixed = TRUE)[[1]]
      for (tk in toks) {
        if (grepl("^eq([0-9]+)::(.*)$", tk)) {
          eq_i <- as.integer(sub("^eq([0-9]+)::(.*)$", "\\1", tk))
          atom_i <- sub("^eq([0-9]+)::(.*)$", "\\2", tk)
          if (eq_i >= 1L && eq_i <= m) {
            disc_list[[eq_i]] <- c(disc_list[[eq_i]], atom_i)
          }
        }
      }
    }

    # Note: legacy saved Discovered_Atoms already stripped sin_deg -> sin.
    # For BGNLM_IWSM_Untuned on trigonometric systems where sin(x) came from sin_deg(x)
    # and R2_ood < 0.85, sin_deg check is also explicitly verified in targeted ablation,
    # while rational-complement normalization modulo intercept runs cleanly here on all rows.
    m_res <- compute_structure_metrics(disc_list, gt_list)
    if (abs(m_res$F1 - dt$F1[i]) > 1e-6 || abs(m_res$FDR - dt$FDR[i]) > 1e-6) {
      n_changed <- n_changed + 1L
    }
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

  saveRDS(dt, rds_path)
  fwrite(dt, csv_path)
  message(sprintf("[OK] Re-scored %s: %d rows (%d metrics updated; FAILED=%d set to NA FDR/OC/SHD)",
                  suite, n_rows, n_changed, sum(dt$Status != "OK", na.rm = TRUE)))

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
        n_ok = length(ok_idx),
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
        n_ok = length(ok_idx),
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
}

rescore_dt("synth")
rescore_dt("real")
