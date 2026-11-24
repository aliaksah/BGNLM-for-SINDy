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
})

source("R/01_ode_simulators.R")
source("R/02_basis_builders.R")
source("R/03_fbms_wrappers.R")
source("R/04_sindy_baselines.R")
source("R/04b_sr_baselines.R")
source("R/05_symbolic_canonicalize.R")
source("R/06_metrics_consistency.R")
source("R/10_real_data_experiments.R")

dir.create("results/ablation_native", recursive = TRUE, showWarnings = FALSE)

# Ensure native SR and fixed standardized UQ-SINDy are active
Sys.setenv(
  PYSR_USE_JULIA = Sys.getenv("PYSR_USE_JULIA", "1"),
  PYOPERON_NATIVE = Sys.getenv("PYOPERON_NATIVE", "1"),
  SINDY_LEGACY_UNSTD_THRESH = "0",
  USE_AUTHENTIC_EMPIRICAL = "1"
)

cores_arg <- as.integer(Sys.getenv("CORES", "4"))
seeds <- c(101L, 202L, 303L)
systems <- c("sys02_lorenz63", "sys08_rossler_lorenz_hybrid", "sys11_goodwin_mm")
noise_regimes <- c("derivative_snr", "state_noise")
methods_synth <- c(
  "BGNLM_Full_Tuned",
  "SINDy_STLSQ_Oracle",
  "SINDy_STLSQ_Poly",
  "UQ_SINDy_Poly",
  "SR_PySR",
  "SR_Operon"
)

run_method_dispatch <- function(method_id, task_data, cores = 4L) {
  if (method_id %in% c("BGNLM_Full_Tuned", "BGLM_Oracle_Basis", "BGLM_Poly_Basis")) {
    run_fbms_method(method_id, task_data, cores = cores)
  } else if (method_id %in% c("SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly", "UQ_SINDy_Oracle", "UQ_SINDy_Poly")) {
    run_sindy_method(method_id, task_data)
  } else {
    run_sr_baseline(method_id, task_data)
  }
}

synth_rows <- list()
for (sys_id in systems) {
  for (reg in noise_regimes) {
    for (sd in seeds) {
      task_data <- tryCatch(
        generate_synthetic_task_data(
          system_id = sys_id,
          N_train = 250L,
          noise_eta = 0.05,
          noise_regime = reg,
          seed = sd
        ),
        error = function(e) {
          message("[SKIP SIM] ", sys_id, " ", reg, ": ", e$message)
          NULL
        }
      )
      if (is.null(task_data)) next
      for (m_id in methods_synth) {
        t0 <- proc.time()[3]
        fit_res <- tryCatch(run_method_dispatch(m_id, task_data, cores = cores_arg), error = function(e) NULL)
        elapsed <- as.numeric(proc.time()[3] - t0)
        if (is.null(fit_res)) {
          synth_rows[[length(synth_rows) + 1L]] <- data.table(
            System = sys_id, NoiseRegime = reg, Seed = sd, Method = m_id,
            Status = "FAILED", F1 = 0, TPR = 0, FDR = 0, Exact = 0,
            R2_ood = 0, R2_traj = 0, Runtime_s = elapsed
          )
        } else {
          ev <- evaluate_method_predictions(fit_res, task_data)
          synth_rows[[length(synth_rows) + 1L]] <- data.table(
            System = sys_id, NoiseRegime = reg, Seed = sd, Method = m_id,
            Status = "OK", F1 = ev$F1, TPR = ev$TPR, FDR = ev$FDR, Exact = ev$Exact_Recovery,
            R2_ood = ev$R2_test_ood, R2_traj = ev$R2_traj, Runtime_s = elapsed
          )
        }
      }
    }
  }
}

dt_synth_abl <- rbindlist(synth_rows, fill = TRUE)
fwrite(dt_synth_abl, "results/ablation_native/ablation_synth_results.csv")
message("[OK] Saved results/ablation_native/ablation_synth_results.csv (", nrow(dt_synth_abl), " rows).")

# Authentic 1900-1920 Hudson Bay Company Lynx-Hare empirical benchmark
if (file.exists("scripts/download_real_data.R")) {
  system("Rscript scripts/download_real_data.R")
}
auth_csv <- "data/real_world/lynx_hare_authentic_1900_1920.csv"
if (file.exists(auth_csv)) {
  raw_df <- fread(auth_csv)
  times <- as.numeric(raw_df$time)
  X_clean <- as.matrix(raw_df[, .(x1, x2)])
  Y_sg <- estimate_derivatives_sgolay_spline(X_clean, times)
  n_tot <- nrow(X_clean)
  n_tr <- floor(0.7 * n_tot)
  tr_idx <- seq_len(n_tr)
  te_idx <- seq(n_tr + 1L, n_tot)
  auth_task <- list(
    system_id = "authentic_lynx_hare_1900_1920",
    X_train = X_clean[tr_idx, , drop = FALSE],
    Y_train = Y_sg[tr_idx, , drop = FALSE],
    X_test_in = X_clean[tr_idx, , drop = FALSE],
    Y_test_in = Y_sg[tr_idx, , drop = FALSE],
    X_test_ood = X_clean[te_idx, , drop = FALSE],
    Y_test_ood = Y_sg[te_idx, , drop = FALSE],
    times_ood = times[te_idx],
    x0_ood = X_clean[te_idx[1L], ],
    ground_truth = list(c("x1", "x1*x2"), c("x2", "x1*x2")),
    seed = 1900L
  )
  real_rows <- list()
  for (m_id in methods_synth) {
    t0 <- proc.time()[3]
    fit_res <- tryCatch(run_method_dispatch(m_id, auth_task, cores = cores_arg), error = function(e) NULL)
    elapsed <- as.numeric(proc.time()[3] - t0)
    if (is.null(fit_res)) {
      real_rows[[length(real_rows) + 1L]] <- data.table(
        Dataset = "authentic_lynx_hare", Method = m_id, Status = "FAILED",
        F1 = 0, FDR = 0, R2_ood = 0, R2_traj = 0, Runtime_s = elapsed
      )
    } else {
      ev <- evaluate_method_predictions(fit_res, auth_task)
      real_rows[[length(real_rows) + 1L]] <- data.table(
        Dataset = "authentic_lynx_hare", Method = m_id, Status = "OK",
        F1 = ev$F1, FDR = ev$FDR, R2_ood = ev$R2_test_ood, R2_traj = ev$R2_traj, Runtime_s = elapsed
      )
    }
  }
  dt_real_abl <- rbindlist(real_rows, fill = TRUE)
  fwrite(dt_real_abl, "results/ablation_native/ablation_authentic_lynx_hare.csv")
  message("[OK] Saved results/ablation_native/ablation_authentic_lynx_hare.csv.")
}
