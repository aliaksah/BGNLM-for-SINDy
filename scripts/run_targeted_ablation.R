#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_targeted_ablation.R
# Targeted Follow-Up Ablation:
#   1. Noise-regime comparison: derivative_snr vs. state_noise (errors-in-variables)
#   2. Native Julia/C++ Symbolic Regression (PYSR_USE_JULIA=1, PYOPERON_NATIVE=1)
#   3. Standardized-threshold UQ_SINDy fix (SINDY_LEGACY_UNSTD_THRESH=0)
#   4. Authentic 1900-1920 Hudson Bay Company Lynx-Hare pelt data
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

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

dir.create("results/ablation_native", recursive = TRUE, showWarnings = FALSE)

# Ensure native SR, fixed standardized UQ-SINDy, and authentic Lynx-Hare are active
Sys.setenv(
  PYSR_USE_JULIA = Sys.getenv("PYSR_USE_JULIA", "1"),
  PYOPERON_NATIVE = Sys.getenv("PYOPERON_NATIVE", "1"),
  SINDY_LEGACY_UNSTD_THRESH = "0",
  USE_AUTHENTIC_EMPIRICAL = "1"
)

smoke_mode <- (Sys.getenv("ABLATION_SMOKE", "0") == "1")
cores_arg <- as.integer(Sys.getenv("CORES", "4"))
timeout_sec <- as.numeric(Sys.getenv("TASK_TIMEOUT", "600"))

seeds <- if (smoke_mode) c(101L) else c(101L, 202L, 303L)
systems <- if (smoke_mode) c("sys02_lorenz63") else c("sys02_lorenz63", "sys08_rossler_lorenz_hybrid", "sys11_goodwin_mm")
noise_regimes <- c("derivative_snr", "state_noise")
methods_synth <- c(
  "BGNLM_Full_Tuned",
  "SINDy_STLSQ_Oracle",
  "SINDy_STLSQ_Poly",
  "UQ_SINDy_Poly",
  "SR_PySR",
  "SR_Operon"
)

run_method_dispatch <- function(method_id, task_data, cores = 4L, fast_mode = FALSE) {
  if (method_id %in% c(
    "BGNLM_Full_Tuned", "BGNLM_IWSM_Untuned", "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis", "BGLM_Poly_Basis", "BGLM_Overcomplete"
  )) {
    run_fbms_method(method_id = method_id, task_data = task_data, cores = cores, fast_mode = fast_mode)
  } else if (method_id %in% c(
    "SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly",
    "UQ_SINDy_Oracle", "UQ_SINDy_Poly", "SINDy_PI_Rational"
  )) {
    run_sindy_baseline(method_id = method_id, task_data = task_data, n_boot = if (fast_mode) 20L else 80L)
  } else if (method_id %in% c("SR_PySR", "SR_Operon")) {
    run_sr_baseline(method_id = method_id, task_data = task_data)
  } else {
    stop("Unknown method_id: ", method_id)
  }
}

score_method_fit <- function(fit_res, task_data, method_id) {
  if (is.null(fit_res)) {
    return(list(
      Status = "FAILED", Engine = "failed", F1 = 0.0, TPR = 0.0, FDR = 0.0, Exact = 0L,
      SHD = NA_integer_, OC = 0L, PIP_Brier = NA_real_,
      R2_ood = 0.0, R2_traj = 0.0, Coverage_95_ood = NA_real_
    ))
  }
  eng <- if (!is.null(fit_res$engine)) as.character(fit_res$engine)[1L] else "r_native"
  st_m <- compute_structure_metrics(fit_res$discovered, task_data$system$ground_truth)
  pip_b <- compute_pip_brier_score(fit_res$pip_list, task_data$system$ground_truth)
  pred_ood_noisy <- compute_prediction_metrics(fit_res$pred_ood, task_data$Y_test_ood_noisy, fit_res$lower_ood, fit_res$upper_ood)
  pred_ood_clean <- compute_prediction_metrics(fit_res$pred_ood, task_data$Y_test_ood_clean)
  roll_m <- evaluate_trajectory_rollout(fit_res$predict_rhs_fn, task_data$rollout)
  list(
    Status = "OK",
    Engine = eng,
    F1 = st_m$F1,
    TPR = st_m$TPR,
    FDR = st_m$FDR,
    Exact = st_m$Exact_Recovery,
    SHD = st_m$SHD,
    OC = st_m$OC,
    PIP_Brier = pip_b,
    R2_ood = pred_ood_clean$R2_clipped,
    R2_traj = roll_m$R2_traj,
    Coverage_95_ood = pred_ood_noisy$Coverage_95
  )
}

build_authentic_lynx_hare_task <- function(csv_path, diff_method = "sgolay_spline", fold_id = 1L, seed = 1900L) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  times <- as.numeric(df$time)
  state_names <- c("x1", "x2")
  X_states <- as.matrix(df[, state_names, drop = FALSE])
  Y_est <- if (diff_method == "sgolay_spline") {
    estimate_derivatives_sgolay_spline(X_states, times)
  } else {
    estimate_derivatives_tvregdiff(X_states, times)
  }
  n_total <- nrow(X_states)
  set.seed(as.integer(seed))
  frac_train <- 0.65 + 0.10 * ((fold_id - 1L) %% 10L) / 9.0
  n_train_window <- max(20L, floor(frac_train * n_total))
  idx_train_pool <- seq_len(n_train_window)
  idx_ood <- seq(from = n_train_window + 1L, to = n_total)
  n_sub <- max(15L, floor(0.85 * length(idx_train_pool)))
  idx_train <- sort(sample(idx_train_pool, size = n_sub, replace = FALSE))
  idx_in <- setdiff(idx_train_pool, idx_train)
  if (length(idx_in) < 3L) idx_in <- idx_train
  rollout_idx <- idx_ood[seq_len(min(length(idx_ood), 100L))]
  rollout_times <- times[rollout_idx] - times[rollout_idx[1]]

  list(
    system = list(
      sys_id = 101L,
      id = "authentic_lynx_hare_1900_1920",
      name = "Authentic 1900-1920 Hudson Bay Lynx-Hare Pelts",
      group = "Ecology",
      dim = 2L,
      state_names = state_names,
      feature_names = state_names,
      input_names = character(0),
      ground_truth = list(c("x1", "x1*x2"), c("x2", "x1*x2")),
      oracle_extra_basis = character(0)
    ),
    dataset_id = "authentic_lynx_hare_1900_1920",
    diff_method = diff_method,
    fold_id = fold_id,
    seed = seed,
    N_train = length(idx_train),
    X_train = X_states[idx_train, , drop = FALSE],
    Y_train = Y_est[idx_train, , drop = FALSE],
    Y_train_clean = Y_est[idx_train, , drop = FALSE],
    X_test_in = X_states[idx_in, , drop = FALSE],
    Y_test_in_clean = Y_est[idx_in, , drop = FALSE],
    Y_test_in_noisy = Y_est[idx_in, , drop = FALSE],
    X_test_ood = X_states[idx_ood, , drop = FALSE],
    Y_test_ood_clean = Y_est[idx_ood, , drop = FALSE],
    Y_test_ood_noisy = Y_est[idx_ood, , drop = FALSE],
    rollout = list(
      times = rollout_times,
      X_true = X_states[rollout_idx, , drop = FALSE],
      X_inputs = NULL,
      x0 = X_states[rollout_idx[1], ]
    )
  )
}

# 1. Synthetic derivative_snr vs. state_noise ablation
synth_rows <- list()
total_synth <- length(systems) * length(noise_regimes) * length(seeds) * length(methods_synth)
idx_run <- 0L

for (sys_id in systems) {
  for (reg in noise_regimes) {
    for (sd in seeds) {
      task_data <- generate_synthetic_task_data(
        system_id = sys_id,
        N_train = if (smoke_mode) 80L else 250L,
        noise_eta = 0.05,
        noise_regime = reg,
        seed = sd,
        T_max = if (smoke_mode) 10.0 else NULL
      )
      for (m_id in methods_synth) {
        idx_run <- idx_run + 1L
        t0 <- proc.time()[3]
        fit_res <- tryCatch(
          R.utils::withTimeout(
            run_method_dispatch(m_id, task_data, cores = cores_arg, fast_mode = smoke_mode),
            timeout = timeout_sec, cpu = Inf, elapsed = timeout_sec, onTimeout = "error"
          ),
          error = function(e) {
            message("[WARN] ", sys_id, " | ", reg, " | ", m_id, " failed: ", e$message)
            NULL
          }
        )
        elapsed <- as.numeric(proc.time()[3] - t0)
        ev <- score_method_fit(fit_res, task_data, m_id)
        message(sprintf(
          "[SYNTH %d/%d] %s | %s | seed=%d | %s (%s) -> %s (F1=%.3f, FDR=%.3f, R2_traj=%.3f, %.2fs)",
          idx_run, total_synth, sys_id, reg, sd, m_id, ev$Engine, ev$Status, ev$F1, ev$FDR, ev$R2_traj, elapsed
        ))
        synth_rows[[length(synth_rows) + 1L]] <- data.table(
          System = sys_id, NoiseRegime = reg, Seed = sd, Method = m_id, Engine = ev$Engine,
          Status = ev$Status, F1 = ev$F1, TPR = ev$TPR, FDR = ev$FDR, Exact = ev$Exact,
          OC = ev$OC, PIP_Brier = ev$PIP_Brier,
          R2_ood = ev$R2_ood, R2_traj = ev$R2_traj, Coverage_95_ood = ev$Coverage_95_ood,
          Runtime_s = elapsed
        )
      }
    }
  }
}

dt_synth_abl <- rbindlist(synth_rows, fill = TRUE)
out_synth_csv <- if (smoke_mode) "results/ablation_native/smoke_synth.csv" else "results/ablation_native/ablation_synth_results.csv"
fwrite(dt_synth_abl, out_synth_csv)
message("[OK] Saved ", out_synth_csv, " (", nrow(dt_synth_abl), " rows).")

# 2. Authentic 1900-1920 Hudson Bay Company Lynx-Hare empirical benchmark
source("scripts/download_real_data.R")
auth_csv <- "data/real_world/lynx_hare_authentic_1900_1920.csv"
if (file.exists(auth_csv)) {
  folds <- if (smoke_mode) 1L else 1:5
  real_rows <- list()
  for (fid in folds) {
    auth_task <- build_authentic_lynx_hare_task(auth_csv, diff_method = "sgolay_spline", fold_id = fid, seed = 1900L + fid)
    for (m_id in methods_synth) {
      t0 <- proc.time()[3]
      fit_res <- tryCatch(
        R.utils::withTimeout(
          run_method_dispatch(m_id, auth_task, cores = cores_arg, fast_mode = smoke_mode),
          timeout = timeout_sec, cpu = Inf, elapsed = timeout_sec, onTimeout = "error"
        ),
        error = function(e) {
          message("[WARN] authentic_lynx_hare | fold=", fid, " | ", m_id, " failed: ", e$message)
          NULL
        }
      )
      elapsed <- as.numeric(proc.time()[3] - t0)
      ev <- score_method_fit(fit_res, auth_task, m_id)
      message(sprintf(
        "[AUTH REAL fold=%d] %s (%s) -> %s (F1=%.3f, FDR=%.3f, R2_traj=%.3f, %.2fs)",
        fid, m_id, ev$Engine, ev$Status, ev$F1, ev$FDR, ev$R2_traj, elapsed
      ))
      real_rows[[length(real_rows) + 1L]] <- data.table(
        Dataset = "authentic_lynx_hare_1900_1920", Fold_ID = fid, Method = m_id, Engine = ev$Engine,
        Status = ev$Status, F1 = ev$F1, TPR = ev$TPR, FDR = ev$FDR, Exact = ev$Exact,
        OC = ev$OC, R2_ood = ev$R2_ood, R2_traj = ev$R2_traj, Coverage_95_ood = ev$Coverage_95_ood,
        Runtime_s = elapsed
      )
    }
  }
  dt_real_abl <- rbindlist(real_rows, fill = TRUE)
  out_real_csv <- if (smoke_mode) "results/ablation_native/smoke_real.csv" else "results/ablation_native/ablation_authentic_lynx_hare.csv"
  fwrite(dt_real_abl, out_real_csv)
  message("[OK] Saved ", out_real_csv, " (", nrow(dt_real_abl), " rows).")
}
