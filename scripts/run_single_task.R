#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_single_task.R
# Atomic CLI Worker for Synthetic & Real-World Benchmark Tasks with Checkpointing
# Usage: Rscript scripts/run_single_task.R --suite <synth|real> --task_id <int> [--cores <int>]
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
})

source("R/01_ode_simulators.R")
source("R/02_basis_builders.R")
source("R/05_symbolic_canonicalize.R")
source("R/06_metrics_consistency.R")
source("R/03_fbms_wrappers.R")
source("R/04_sindy_baselines.R")
source("R/04b_sr_baselines.R")
source("R/10_real_data_experiments.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0L && length(args) > idx[1]) return(args[idx[1] + 1L])
  default
}

suite <- get_arg("--suite", "synth")
task_id <- as.integer(get_arg("--task_id", "1"))
cores <- as.integer(get_arg("--cores", "2"))
timeout_sec <- as.numeric(get_arg("--timeout", "600"))

manifest_path <- paste0("config/manifest_", suite, ".csv")
if (!file.exists(manifest_path)) {
  stop("Manifest file does not exist: ", manifest_path)
}

manifest <- fread(manifest_path)
if (task_id < 1L || task_id > nrow(manifest)) {
  stop("Invalid task_id ", task_id, " (manifest has ", nrow(manifest), " rows).")
}

ckpt_dir <- paste0("results/checkpoints_", suite)
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)
ckpt_file <- file.path(ckpt_dir, sprintf("task_%06d.rds", task_id))

# Seamless resumption: if checkpoint already exists, exit cleanly (0)
if (file.exists(ckpt_file)) {
  message("[SKIP] Checkpoint already exists: ", ckpt_file)
  quit(save = "no", status = 0L)
}

row <- manifest[task_id]
method_id <- as.character(row$method_id)
fast_mode <- as.logical(row$fast_mode == 1L)

# Generate deterministic dataset split
if (suite == "synth") {
  T_sim <- if (fast_mode) 15.0 else NULL
  task_data <- generate_synthetic_task_data(
    system_id = as.character(row$system_id),
    N_train = as.integer(row$sample_size),
    noise_eta = as.numeric(row$noise_eta),
    noise_regime = as.character(row$noise_regime),
    seed = as.integer(row$seed),
    T_max = T_sim
  )
} else {
  task_data <- generate_real_task_data(
    dataset_id = as.character(row$dataset_id),
    diff_method = as.character(row$diff_method),
    fold_id = as.integer(row$fold_id),
    seed = as.integer(row$seed)
  )
}

# Execute method under CPU/elapsed timeout guard
t_start <- proc.time()[3]
gc_before <- gc(reset = TRUE)

fit_res <- tryCatch(
  R.utils::withTimeout({
    if (method_id %in% c(
      "BGNLM_Full_Tuned", "BGNLM_IWSM_Untuned", "BGNLM_Hybrid_Poly",
      "BGLM_Oracle_Basis", "BGLM_Poly_Basis", "BGLM_Overcomplete"
    )) {
      run_fbms_method(method_id = method_id, task_data = task_data, cores = cores, fast_mode = fast_mode)
    } else if (method_id %in% c(
      "SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly",
      "UQ_SINDy_Oracle", "UQ_SINDy_Poly", "SINDy_PI_Rational"
    )) {
      run_sindy_baseline(method_id = method_id, task_data = task_data, n_boot = if (fast_mode) 40L else 80L)
    } else if (method_id %in% c("SR_PySR", "SR_Operon")) {
      run_sr_baseline(method_id = method_id, task_data = task_data)
    } else {
      stop("Unknown method_id: ", method_id)
    }
  }, timeout = timeout_sec, onTimeout = "error"),
  error = function(e) {
    warning("Method ", method_id, " encountered error on task ", task_id, ": ", e$message)
    NULL
  }
)

runtime_sec <- as.numeric(proc.time()[3] - t_start)
gc_after <- gc()
peak_mem_mb <- sum(gc_after[, 6], na.rm = TRUE)

# Compute all structural, calibration, predictive, and trajectory rollout metrics
m <- ncol(task_data$Y_train)
if (is.null(fit_res)) {
  discovered_empty <- replicate(m, character(0), simplify = FALSE)
  st_metrics <- compute_structure_metrics(discovered_empty, task_data$system$ground_truth)
  pip_brier <- NA_real_
  pred_tr <- list(R2_clipped = 0.0, R2_raw = 0.0, NRMSE = 1.0, log10_RMSE = 0.0)
  pred_in <- list(R2_clipped = 0.0, R2_raw = 0.0, NRMSE = 1.0, log10_RMSE = 0.0, Coverage_95 = NA_real_, Winkler_95 = NA_real_)
  pred_ood <- list(R2_clipped = 0.0, R2_raw = 0.0, NRMSE = 1.0, log10_RMSE = 0.0, Coverage_95 = NA_real_, Winkler_95 = NA_real_)
  rollout_m <- list(R2_traj = 0.0, NRMSE_traj = 1.0, T_valid = 0.0, T_valid_frac = 0.0)
  status_str <- "FAILED"
} else {
  st_metrics <- compute_structure_metrics(fit_res$discovered, task_data$system$ground_truth)
  pip_brier <- compute_pip_brier_score(fit_res$pip_list, task_data$system$ground_truth)
  pred_tr <- compute_prediction_metrics(fit_res$pred_train, task_data$Y_train_clean)
  pred_in <- compute_prediction_metrics(fit_res$pred_in, task_data$Y_test_in_noisy, fit_res$lower_in, fit_res$upper_in)
  pred_ood <- compute_prediction_metrics(fit_res$pred_ood, task_data$Y_test_ood_noisy, fit_res$lower_ood, fit_res$upper_ood)
  pred_in_clean <- compute_prediction_metrics(fit_res$pred_in, task_data$Y_test_in_clean)
  pred_ood_clean <- compute_prediction_metrics(fit_res$pred_ood, task_data$Y_test_ood_clean)
  pred_in$R2_clipped <- pred_in_clean$R2_clipped
  pred_in$R2_raw <- pred_in_clean$R2_raw
  pred_ood$R2_clipped <- pred_ood_clean$R2_clipped
  pred_ood$R2_raw <- pred_ood_clean$R2_raw
  rollout_m <- evaluate_trajectory_rollout(fit_res$predict_rhs_fn, task_data$rollout)
  status_str <- "OK"
}

if (suite == "synth") {
  out_dt <- data.table(
    task_id = task_id,
    suite = "synth",
    System = as.character(row$system_id),
    SystemGroup = as.character(row$system_group),
    SampleSize = as.integer(row$sample_size),
    Noise = as.numeric(row$noise_eta),
    NoiseRegime = as.character(row$noise_regime),
    Rep_ID = as.integer(row$rep_id),
    Seed = as.integer(row$seed),
    Method = method_id,
    Status = status_str,
    TPR = st_metrics$TPR,
    FDR = st_metrics$FDR,
    F1 = st_metrics$F1,
    Exact_Recovery = st_metrics$Exact_Recovery,
    SHD = st_metrics$SHD,
    OC = st_metrics$OC,
    n_active_terms = st_metrics$n_active_terms,
    PIP_Brier = pip_brier,
    R2_train = pred_tr$R2_clipped,
    R2_test_in = pred_in$R2_clipped,
    R2_test_ood = pred_ood$R2_clipped,
    R2_test_ood_raw = pred_ood$R2_raw,
    NRMSE_ood = pred_ood$NRMSE,
    log10_RMSE_ood = pred_ood$log10_RMSE,
    Coverage_95_in = pred_in$Coverage_95,
    Winkler_95_in = pred_in$Winkler_95,
    Coverage_95_ood = pred_ood$Coverage_95,
    Winkler_95_ood = pred_ood$Winkler_95,
    R2_traj = rollout_m$R2_traj,
    NRMSE_traj = rollout_m$NRMSE_traj,
    T_valid = rollout_m$T_valid,
    T_valid_frac = rollout_m$T_valid_frac,
    Runtime_s = runtime_sec,
    Peak_Mem_MB = peak_mem_mb,
    Discovered_Atoms = paste(st_metrics$S_hat, collapse = ";"),
    True_Atoms = paste(st_metrics$S_star, collapse = ";")
  )
} else {
  out_dt <- data.table(
    task_id = task_id,
    suite = "real",
    Dataset = as.character(row$dataset_id),
    Domain = as.character(row$domain),
    DiffMethod = as.character(row$diff_method),
    Fold_ID = as.integer(row$fold_id),
    Seed = as.integer(row$seed),
    Method = method_id,
    Status = status_str,
    TPR = st_metrics$TPR,
    FDR = st_metrics$FDR,
    F1 = st_metrics$F1,
    Exact_Recovery = st_metrics$Exact_Recovery,
    SHD = st_metrics$SHD,
    OC = st_metrics$OC,
    n_active_terms = st_metrics$n_active_terms,
    PIP_Brier = pip_brier,
    R2_train = pred_tr$R2_clipped,
    R2_test_in = pred_in$R2_clipped,
    R2_test_ood = pred_ood$R2_clipped,
    NRMSE_deriv_ood = pred_ood$NRMSE,
    Coverage_95_ood = pred_ood$Coverage_95,
    Winkler_95_ood = pred_ood$Winkler_95,
    R2_traj = rollout_m$R2_traj,
    NRMSE_traj = rollout_m$NRMSE_traj,
    T_valid = rollout_m$T_valid,
    T_valid_frac = rollout_m$T_valid_frac,
    Runtime_s = runtime_sec,
    Peak_Mem_MB = peak_mem_mb,
    Discovered_Atoms = paste(st_metrics$S_hat, collapse = ";"),
    Reference_Atoms = paste(st_metrics$S_star, collapse = ";")
  )
}

saveRDS(out_dt, ckpt_file)
message(sprintf("[DONE] suite=%s task_id=%d method=%s F1=%.3f R2_ood=%.3f R2_traj=%.3f (%.2fs)",
                suite, task_id, method_id, out_dt$F1, out_dt$R2_test_ood, out_dt$R2_traj, runtime_sec))

