#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_tuning.R
# CLI Entrypoint to Execute Phase 4 Systematic BGNLM Prior & GMJMCMC Tuning
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")
source("R/07_tuning_engine.R")

args <- commandArgs(trailingOnly = TRUE)
full_mode <- "--full" %in% args
cores_idx <- which(args == "--cores")
cores <- if (length(cores_idx) > 0L && length(args) > cores_idx) as.integer(args[cores_idx + 1L]) else 2L

res <- execute_bgnlm_tuning(cores = cores, fast_calibration = !full_mode)
print(res$summary[, .(config_id, beta_prior_type, r_rule, p_proj, Mean_F1, Mean_FDR, Exact_Rate, Mean_PIP_Brier, Mean_R2_ood, Composite_Score)])

