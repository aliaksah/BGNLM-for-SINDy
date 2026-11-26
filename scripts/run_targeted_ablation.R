#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_targeted_ablation.R
# Targeted Follow-Up Ablation & Fair Native Benchmark Suite (Jobs 1-4):
#   Job 1: Pure Bayesian MPM (tau=0) vs Pruned MPM (tau=0.025), BMA vs MPM-OLS,
#          Symmetric SINDy pruning, & Large-N (N=200..2500) scaling diagnosis/fix
#   Job 2: Fair Baselines Head-to-Head (Fixed Bag-STLSQ with lambda*, UQ_SINDy_Bayes
#          via BAS::bas.lm, WSINDy Galerkin integral, SINDy_PI_Rat deg-3,
#          Native PySR & PyOperon, sin_deg check on IWSM, & PIP calibration bins)
#   Job 3: Non-Unit Internal Parameter Expressivity Stress Test (sin(2*x1), Km=0.3)
#   Job 4: Zero-Leakage Authentic Real-World Empirical Benchmarks
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

Sys.setenv(
  PYSR_USE_JULIA = Sys.getenv("PYSR_USE_JULIA", "1"),
  PYOPERON_NATIVE = Sys.getenv("PYOPERON_NATIVE", "1"),
  SINDY_LEGACY_UNSTD_THRESH = "0",
  USE_AUTHENTIC_EMPIRICAL = "1"
)

smoke_mode <- (Sys.getenv("ABLATION_SMOKE", "0") == "1")
cores_arg <- as.integer(Sys.getenv("CORES", "2"))
timeout_sec <- as.numeric(Sys.getenv("TASK_TIMEOUT", "420"))
jobs_run <- Sys.getenv("ABLATION_JOBS", "1,2,3,4")
has_job <- function(j_id) grepl(j_id, jobs_run, fixed = TRUE)

run_safe_call <- function(expr_fn) {
  tryCatch(
    R.utils::withTimeout(expr_fn(), timeout = timeout_sec, cpu = Inf, elapsed = timeout_sec, onTimeout = "error"),
    error = function(e) {
      message("[WARN] Fit failed: ", e$message)
      NULL
    }
  )
}

score_fit_full <- function(fit_res, task_data, method_id, discovered_override = NULL, pred_ood_override = NULL) {
  if (is.null(fit_res)) {
    return(list(
      Status = "FAILED", Engine = "failed",
      F1 = 0.0, TPR = 0.0, FDR = NA_real_, Exact = 0L,
      SHD = NA_integer_, OC = NA_integer_,
      PIP_Brier = NA_real_, PIP_SSPE = NA_real_, PIP_AUROC = NA_real_, PIP_AUPRC = NA_real_,
      R2_ood = 0.0, R2_traj = 0.0, Coverage_95_ood = NA_real_,
      Discovered_Atoms = ""
    ))
  }
  eng <- if (!is.null(fit_res$engine)) as.character(fit_res$engine)[1L] else "r_native"
  disc_use <- if (!is.null(discovered_override)) discovered_override else fit_res$discovered
  st_m <- compute_structure_metrics(disc_use, task_data$system$ground_truth)
  cal_m <- compute_pip_calibration_metrics(fit_res$pip_list, task_data$system$ground_truth)
  po_mat <- if (!is.null(pred_ood_override)) pred_ood_override else fit_res$pred_ood
  pred_ood_noisy <- compute_prediction_metrics(po_mat, task_data$Y_test_ood_noisy, fit_res$lower_ood, fit_res$upper_ood)
  pred_ood_clean <- compute_prediction_metrics(po_mat, task_data$Y_test_ood_clean)
  roll_m <- if (is.null(discovered_override) && is.null(pred_ood_override)) {
    evaluate_trajectory_rollout(fit_res$predict_rhs_fn, task_data$rollout)
  } else {
    list(R2_traj = NA_real_)
  }
  list(
    Status = "OK", Engine = eng,
    F1 = st_m$F1, TPR = st_m$TPR, FDR = st_m$FDR, Exact = st_m$Exact_Recovery,
    SHD = st_m$SHD, OC = st_m$OC,
    PIP_Brier = cal_m$PIP_Brier, PIP_SSPE = cal_m$PIP_SSPE,
    PIP_AUROC = cal_m$PIP_AUROC, PIP_AUPRC = cal_m$PIP_AUPRC,
    R2_ood = pred_ood_clean$R2_clipped, R2_traj = roll_m$R2_traj,
    Coverage_95_ood = pred_ood_noisy$Coverage_95,
    Discovered_Atoms = paste(st_m$S_hat, collapse = ";"),
    cal_m = cal_m
  )
}

if (has_job("1")) {
message("=== [JOB 1] MPM Energy Pruning, Symmetric SINDy Pruning, & Large-N Scaling ===")
job1_rows <- list()

# 1a: Pruning & BMA comparison at N=400, sigma=0.05
j1a_systems <- if (smoke_mode) "sys02_lorenz63" else c("sys02_lorenz63", "sys08_rossler_lorenz_hybrid", "sys11_goodwin_mm")
j1a_seeds <- if (smoke_mode) 101L else c(101L, 202L, 303L)

for (sys_id in j1a_systems) {
  for (sd in j1a_seeds) {
    td <- generate_synthetic_task_data(sys_id, N_train = 400L, noise_eta = 0.05, noise_regime = "derivative_snr", seed = sd)

    # Fit BGNLM_Full_Tuned once; extracts BOTH pruned MPM (tau=0.025), Pure MPM (tau=0), and BMA mean
    t0 <- proc.time()[3]
    fit_b <- run_safe_call(function() run_fbms_method("BGNLM_Full_Tuned", td, cores = cores_arg, tuning_override = list(runs = 4L), fast_mode = smoke_mode))
    el_b <- as.numeric(proc.time()[3] - t0)
    ev_pruned <- score_fit_full(fit_b, td, "BGNLM_Full_Tuned")
    ev_pure <- score_fit_full(fit_b, td, "BGNLM_Full_PureMPM", discovered_override = if (!is.null(fit_b)) fit_b$discovered_raw else NULL)
    ev_bma <- if (!is.null(fit_b) && all(is.finite(fit_b$pred_bma_ood))) {
      compute_prediction_metrics(fit_b$pred_bma_ood, td$Y_test_ood_clean)$R2_clipped
    } else NA_real_

    job1_rows[[length(job1_rows) + 1L]] <- data.table(
      SubJob = "1a_Pruning", System = sys_id, SampleSize = 400L, Seed = sd,
      Variant = "BGNLM_Full_Pruned_0.025", F1 = ev_pruned$F1, FDR = ev_pruned$FDR, Exact = ev_pruned$Exact,
      R2_ood = ev_pruned$R2_ood, R2_ood_BMA = ev_bma, R2_traj = ev_pruned$R2_traj, Runtime_s = el_b
    )
    job1_rows[[length(job1_rows) + 1L]] <- data.table(
      SubJob = "1a_Pruning", System = sys_id, SampleSize = 400L, Seed = sd,
      Variant = "BGNLM_Full_PureMPM_0", F1 = ev_pure$F1, FDR = ev_pure$FDR, Exact = ev_pure$Exact,
      R2_ood = ev_pruned$R2_ood, R2_ood_BMA = ev_bma, R2_traj = NA_real_, Runtime_s = el_b
    )

    # Symmetric SINDy pruning comparison (tau=0 vs tau=0.025)
    for (s_mid in c("SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly")) {
      for (ptol in c(0.0, 0.025)) {
        t0s <- proc.time()[3]
        fs <- run_safe_call(function() run_sindy_baseline(s_mid, td, prune_tol = ptol))
        els <- as.numeric(proc.time()[3] - t0s)
        evs <- score_fit_full(fs, td, s_mid)
        v_nm <- paste0(s_mid, ifelse(ptol > 0, "_Pruned_0.025", "_Unpruned"))
        job1_rows[[length(job1_rows) + 1L]] <- data.table(
          SubJob = "1a_Pruning", System = sys_id, SampleSize = 400L, Seed = sd,
          Variant = v_nm, F1 = evs$F1, FDR = evs$FDR, Exact = evs$Exact,
          R2_ood = evs$R2_ood, R2_ood_BMA = NA_real_, R2_traj = evs$R2_traj, Runtime_s = els
        )
      }
    }
    message(sprintf("[JOB 1a] %s seed=%d | Pruned F1=%.3f PureMPM F1=%.3f (BMA R2=%.3f vs MPM-OLS R2=%.3f)",
                    sys_id, sd, ev_pruned$F1, ev_pure$F1, ifelse(is.na(ev_bma), -1, ev_bma), ev_pruned$R2_ood))
  }
}

# 1b: Large-N scaling diagnosis & fix (Fixed budget r=1/N, P=15, N=500 vs Scaled budget r=1/sqrt(N), P=16, N=800)
j1b_specs <- if (smoke_mode) {
  list(list(sys = "sys11_goodwin_mm", N = 200L), list(sys = "sys11_goodwin_mm", N = 1000L))
} else {
  list(
    list(sys = "sys11_goodwin_mm", N = 500L),
    list(sys = "sys11_goodwin_mm", N = 1000L),
    list(sys = "sys11_goodwin_mm", N = 2500L),
    list(sys = "sys02_lorenz63",   N = 2500L)
  )
}
j1b_seeds <- if (smoke_mode) 101L else c(101L, 202L)

for (sp in j1b_specs) {
  sys_id <- sp$sys; N_tr <- sp$N
  for (sd in j1b_seeds) {
    td <- generate_synthetic_task_data(sys_id, N_train = N_tr, noise_eta = 0.05, noise_regime = "derivative_snr", seed = sd)

    # Baseline tuned budget
    t0 <- proc.time()[3]
    f_base <- run_safe_call(function() run_fbms_method("BGNLM_Full_Tuned", td, cores = cores_arg, tuning_override = list(runs = 4L), fast_mode = smoke_mode))
    ev_base <- score_fit_full(f_base, td, "BGNLM_FixedBudget")
    job1_rows[[length(job1_rows) + 1L]] <- data.table(
      SubJob = "1b_LargeN", System = sys_id, SampleSize = N_tr, Seed = sd,
      Variant = "BGNLM_Baseline_FixedBudget", F1 = ev_base$F1, FDR = ev_base$FDR, Exact = ev_base$Exact,
      R2_ood = ev_base$R2_ood, R2_ood_BMA = NA_real_, R2_traj = ev_base$R2_traj, Runtime_s = as.numeric(proc.time()[3] - t0)
    )

    # Scaled budget + sqrt(N) prior fix
    t1 <- proc.time()[3]
    ovr <- list(P_generations = if (smoke_mode) 8L else 16L, N_mcmc = if (smoke_mode) 300L else 800L, runs = 4L, model_prior = list(r = 1.0 / sqrt(N_tr)))
    f_scaled <- run_safe_call(function() run_fbms_method("BGNLM_Full_Tuned", td, cores = cores_arg, tuning_override = ovr, fast_mode = FALSE))
    ev_sc <- score_fit_full(f_scaled, td, "BGNLM_ScaledBudget")
    job1_rows[[length(job1_rows) + 1L]] <- data.table(
      SubJob = "1b_LargeN", System = sys_id, SampleSize = N_tr, Seed = sd,
      Variant = "BGNLM_ScaledBudget_SqrtN", F1 = ev_sc$F1, FDR = ev_sc$FDR, Exact = ev_sc$Exact,
      R2_ood = ev_sc$R2_ood, R2_ood_BMA = NA_real_, R2_traj = ev_sc$R2_traj, Runtime_s = as.numeric(proc.time()[3] - t1)
    )
    message(sprintf("[JOB 1b] %s N=%d seed=%d | Base F1=%.3f | Scaled F1=%.3f", sys_id, N_tr, sd, ev_base$F1, ev_sc$F1))
    fwrite(rbindlist(job1_rows, fill = TRUE), "results/ablation_native/job1_mpm_pruning_and_largeN.csv")
  }
}

dt_job1 <- rbindlist(job1_rows, fill = TRUE)
fwrite(dt_job1, "results/ablation_native/job1_mpm_pruning_and_largeN.csv")
message("[OK] Saved results/ablation_native/job1_mpm_pruning_and_largeN.csv")
}

# ==============================================================================
# JOB 2: Fair Baselines Head-to-Head + sin_deg Check + PIP Calibration Reliability
# ==============================================================================
if (has_job("2")) {
message("=== [JOB 2] Fair Baselines Head-to-Head & PIP Calibration ===")
j2_systems <- if (smoke_mode) c("sys02_lorenz63", "sys09_damped_pendulum") else c("sys02_lorenz63", "sys08_rossler_lorenz_hybrid", "sys09_damped_pendulum", "sys11_goodwin_mm")
j2_regimes <- c("derivative_snr", "state_noise")
j2_seeds <- if (smoke_mode) 101L else c(101L, 202L)
j2_methods <- c(
  "BGNLM_Full_Tuned", "BGNLM_IWSM_Untuned",
  "SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly",
  "UQ_SINDy_Poly", "UQ_SINDy_Bayes_Poly", "UQ_SINDy_Bayes_Oracle",
  "WSINDy_Poly", "WSINDy_Oracle", "SINDy_PI_Rational",
  "SR_PySR", "SR_Operon"
)

job2_rows <- list()
pip_calib_rows <- list()

dispatch_j2 <- function(mid, td) {
  if (mid %in% c("BGNLM_Full_Tuned", "BGNLM_IWSM_Untuned")) {
    run_fbms_method(mid, td, cores = cores_arg, tuning_override = list(runs = 4L), fast_mode = smoke_mode)
  } else if (mid %in% c("SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly", "UQ_SINDy_Poly", "SINDy_PI_Rational")) {
    run_sindy_baseline(mid, td, n_boot = if (smoke_mode) 25L else 80L)
  } else if (mid == "UQ_SINDy_Bayes_Poly") {
    run_uq_sindy_bayes(td, dict_type = "poly")
  } else if (mid == "UQ_SINDy_Bayes_Oracle") {
    run_uq_sindy_bayes(td, dict_type = "oracle")
  } else if (mid == "WSINDy_Poly") {
    run_wsindy_baseline(td, dict_type = "poly")
  } else if (mid == "WSINDy_Oracle") {
    run_wsindy_baseline(td, dict_type = "oracle")
  } else if (mid %in% c("SR_PySR", "SR_Operon")) {
    run_sr_baseline(mid, td, strict_native = TRUE)
  } else stop("Unknown method: ", mid)
}

for (sys_id in j2_systems) {
  base_tr <- simulate_base_trajectory(sys_id)
  for (reg in j2_regimes) {
    for (sd in j2_seeds) {
      td <- generate_synthetic_task_data(sys_id, N_train = 300L, noise_eta = 0.05, noise_regime = reg, seed = sd)
      if (reg == "state_noise") {
        # Provide N=300 equispaced contiguous noisy state trajectory from [0.05T, 0.55T] for WSINDy (same N & SNR, zero test leakage)
        idx_win <- which(base_tr$times > 0.05 * base_tr$T_max & base_tr$times < 0.55 * base_tr$T_max)
        idx_eq <- idx_win[round(seq(1L, length(idx_win), length.out = 300L))]
        X_eq <- base_tr$X[idx_eq, , drop = FALSE]
        set.seed(sd + 999L)
        X_eq_noisy <- X_eq + sapply(seq_len(ncol(X_eq)), function(c_i) rnorm(nrow(X_eq), 0.0, td$noise_sd_vec[c_i]))
        td$t_seq <- base_tr$times[idx_eq]
        td$X_seq <- X_eq_noisy
      }
      for (mid in j2_methods) {
        t0 <- proc.time()[3]
        fr <- run_safe_call(function() dispatch_j2(mid, td))
        el <- as.numeric(proc.time()[3] - t0)
        ev <- score_fit_full(fr, td, mid)

        if (mid == "BGNLM_Full_Tuned" && !is.null(fr)) {
          ev_raw <- score_fit_full(fr, td, "BGNLM_Full_PureMPM", discovered_override = fr$discovered_raw)
          job2_rows[[length(job2_rows) + 1L]] <- data.table(
            System = sys_id, NoiseRegime = reg, Seed = sd, Method = "BGNLM_Full_PureMPM", Engine = ev$Engine,
            Status = "OK", F1 = ev_raw$F1, TPR = ev_raw$TPR, FDR = ev_raw$FDR, Exact = ev_raw$Exact,
            SHD = ev_raw$SHD, OC = ev_raw$OC, PIP_Brier = ev$PIP_Brier, PIP_SSPE = ev$PIP_SSPE,
            PIP_AUROC = ev$PIP_AUROC, PIP_AUPRC = ev$PIP_AUPRC,
            R2_ood = ev$R2_ood, R2_traj = NA_real_, Coverage_95_ood = ev$Coverage_95_ood,
            Runtime_s = el, Discovered_Atoms = ev_raw$Discovered_Atoms
          )
        }

        if (!is.null(ev$cal_m) && nrow(ev$cal_m$reliability) > 0L &&
            mid %in% c("BGNLM_Full_Tuned", "UQ_SINDy_Bayes_Poly", "UQ_SINDy_Poly")) {
          rel_dt <- as.data.table(ev$cal_m$reliability)
          rel_dt[, `:=`(System = sys_id, NoiseRegime = reg, Seed = sd, Method = mid)]
          pip_calib_rows[[length(pip_calib_rows) + 1L]] <- rel_dt
        }

        job2_rows[[length(job2_rows) + 1L]] <- data.table(
          System = sys_id, NoiseRegime = reg, Seed = sd, Method = mid, Engine = ev$Engine,
          Status = ev$Status, F1 = ev$F1, TPR = ev$TPR, FDR = ev$FDR, Exact = ev$Exact,
          SHD = ev$SHD, OC = ev$OC, PIP_Brier = ev$PIP_Brier, PIP_SSPE = ev$PIP_SSPE,
          PIP_AUROC = ev$PIP_AUROC, PIP_AUPRC = ev$PIP_AUPRC,
          R2_ood = ev$R2_ood, R2_traj = ev$R2_traj, Coverage_95_ood = ev$Coverage_95_ood,
          Runtime_s = el, Discovered_Atoms = ev$Discovered_Atoms
        )
        message(sprintf("[JOB 2] %s | %s | seed=%d | %s (%s) -> F1=%.3f FDR=%.3f R2_traj=%.3f (%.1fs)",
                        sys_id, reg, sd, mid, ev$Engine, ev$F1, ifelse(is.na(ev$FDR), -1, ev$FDR), ev$R2_traj, el))
      }
    }
  }
}

fwrite(rbindlist(job2_rows, fill = TRUE), "results/ablation_native/job2_fair_baselines_synth.csv")
if (length(pip_calib_rows) > 0L) {
  fwrite(rbindlist(pip_calib_rows, fill = TRUE), "results/ablation_native/pip_calibration_reliability.csv")
}
message("[OK] Saved job2_fair_baselines_synth.csv & pip_calibration_reliability.csv")
}

# ==============================================================================
# JOB 3: Non-Unit Internal Parameter Expressivity Stress Test
#   System A: Forced Non-Unit Pendulum dx1 = x2, dx2 = -0.25*x2 - 5*sin(2*x1)
#   System B: Non-Unit Km Goodwin dx1 = 1/(1+x3^2) - 0.5*x1, dx2 = x1 - 0.8*x2/(0.3 + x2), dx3 = x2 - 0.5*x3
# ==============================================================================
if (has_job("3")) {
message("=== [JOB 3] Non-Unit Internal Parameter Stress Test ===")
build_custom_task <- function(sys_def, N_train = 300L, noise_eta = 0.03, seed = 101L) {
  times <- seq(0, sys_def$T_max, by = sys_def$dt)
  sol <- as.data.frame(deSolve::ode(y = sys_def$y0, times = times, func = sys_def$rhs, parms = sys_def$parms, method = "lsoda"))
  X_all <- as.matrix(sol[, -1, drop = FALSE]); colnames(X_all) <- sys_def$state_names
  t_all <- sol[, 1]
  dX_all <- sys_def$analytical_deriv(X_all, sys_def$parms)
  set.seed(seed)
  n_tot <- nrow(X_all)
  idx_int <- which(t_all > 0.05 * sys_def$T_max & t_all < 0.90 * sys_def$T_max)
  idx_ood_all <- which(t_all >= 0.90 * sys_def$T_max)
  samp <- sort(sample(idx_int, size = min(length(idx_int), N_train + 300L)))
  idx_tr <- samp[seq_len(N_train)]
  idx_in <- setdiff(samp, idx_tr)
  idx_ood <- idx_ood_all[round(seq(1, length(idx_ood_all), length.out = 300L))]

  m <- ncol(X_all)
  sd_d <- pmax(apply(dX_all[idx_int, , drop = FALSE], 2, sd), 1e-6)
  eps_tr <- sapply(seq_len(m), function(j) rnorm(length(idx_tr), 0, noise_eta * sd_d[j]))
  eps_in <- sapply(seq_len(m), function(j) rnorm(length(idx_in), 0, noise_eta * sd_d[j]))
  eps_ood <- sapply(seq_len(m), function(j) rnorm(length(idx_ood), 0, noise_eta * sd_d[j]))
  roll_idx <- idx_ood_all[seq(1L, min(length(idx_ood_all), 2000L), by = 10L)]

  list(
    system = list(
      sys_id = sys_def$sys_id, id = sys_def$id, name = sys_def$name,
      dim = m, state_names = sys_def$state_names, feature_names = sys_def$state_names,
      input_names = character(0), ground_truth = sys_def$ground_truth,
      oracle_extra_basis = character(0)
    ),
    seed = seed, N_train = N_train,
    t_train = t_all[idx_tr],
    X_train = X_all[idx_tr, , drop = FALSE], Y_train = dX_all[idx_tr, , drop = FALSE] + eps_tr,
    Y_train_clean = dX_all[idx_tr, , drop = FALSE],
    X_test_in = X_all[idx_in, , drop = FALSE],
    Y_test_in_clean = dX_all[idx_in, , drop = FALSE], Y_test_in_noisy = dX_all[idx_in, , drop = FALSE] + eps_in,
    X_test_ood = X_all[idx_ood, , drop = FALSE],
    Y_test_ood_clean = dX_all[idx_ood, , drop = FALSE], Y_test_ood_noisy = dX_all[idx_ood, , drop = FALSE] + eps_ood,
    rollout = list(times = t_all[roll_idx] - t_all[roll_idx[1]], X_true = X_all[roll_idx, , drop = FALSE], x0 = X_all[roll_idx[1], ])
  )
}

nonunit_systems <- list(
  nonunit_pendulum_2theta = list(
    sys_id = 901L, id = "nonunit_pendulum_2theta", name = "Non-Unit Frequency Pendulum sin(2*x1)",
    dt = 0.002, T_max = 20.0, y0 = c(x1 = 1.8, x2 = 0.0), parms = c(g = 0.25, a = 5.0),
    state_names = c("x1", "x2"),
    rhs = function(t, s, p) list(c(s[2], -p["g"] * s[2] - p["a"] * sin(2.0 * s[1]))),
    analytical_deriv = function(X, p) cbind(dx1 = X[, 2], dx2 = -p["g"] * X[, 2] - p["a"] * sin(2.0 * X[, 1])),
    ground_truth = list(c("x2"), c("x2", "sin(2*x1)"))
  ),
  nonunit_mm_km03 = list(
    sys_id = 902L, id = "nonunit_mm_km03", name = "Non-Unit Michaelis-Menten (Km=0.3)",
    dt = 0.002, T_max = 30.0, y0 = c(x1 = 0.5, x2 = 0.5, x3 = 0.5), parms = c(k = 0.5),
    state_names = c("x1", "x2", "x3"),
    rhs = function(t, s, p) list(c(1.0 / (1.0 + s[3]^2) - 0.5 * s[1], s[1] - 0.8 * s[2] / (0.3 + abs(s[2])), s[2] - 0.5 * s[3])),
    analytical_deriv = function(X, p) cbind(dx1 = 1.0 / (1.0 + X[, 3]^2) - 0.5 * X[, 1], dx2 = X[, 1] - 0.8 * X[, 2] / (0.3 + abs(X[, 2])), dx3 = X[, 2] - 0.5 * X[, 3]),
    ground_truth = list(c("x1", "1/(1+x3^2)"), c("x1", "x2/(0.3+x2)"), c("x2", "x3"))
  )
)

j3_seeds <- if (smoke_mode) 101L else c(101L, 202L)
job3_rows <- list()

for (s_nm in names(nonunit_systems)) {
  for (sd in j3_seeds) {
    td <- build_custom_task(nonunit_systems[[s_nm]], N_train = 300L, noise_eta = 0.03, seed = sd)

    variants <- list(
      list(id = "BGNLM_Full_Tuned_pproj0", fn = function() run_fbms_method("BGNLM_Full_Tuned", td, cores = cores_arg, fast_mode = smoke_mode)),
      list(id = "BGNLM_Full_Proj005", fn = function() run_fbms_method("BGNLM_Full_Tuned", td, cores = cores_arg, tuning_override = list(p_inter = 0.50, p_mod = 0.35, p_proj = 0.05, p_new = 0.10), fast_mode = smoke_mode)),
      list(id = "SR_PySR", fn = function() run_sr_baseline("SR_PySR", td, strict_native = TRUE)),
      list(id = "SR_Operon", fn = function() run_sr_baseline("SR_Operon", td, strict_native = TRUE))
    )
    for (v in variants) {
      t0 <- proc.time()[3]
      fr <- run_safe_call(v$fn)
      el <- as.numeric(proc.time()[3] - t0)
      ev <- score_fit_full(fr, td, v$id)
      job3_rows[[length(job3_rows) + 1L]] <- data.table(
        System = s_nm, Seed = sd, Method = v$id, Engine = ev$Engine, Status = ev$Status,
        F1 = ev$F1, FDR = ev$FDR, R2_ood = ev$R2_ood, R2_traj = ev$R2_traj,
        OC = ev$OC, Runtime_s = el, Discovered_Atoms = ev$Discovered_Atoms
      )
      message(sprintf("[JOB 3] %s | seed=%d | %s (%s) -> F1=%.3f R2_ood=%.3f R2_traj=%.3f (%.1fs)",
                      s_nm, sd, v$id, ev$Engine, ev$F1, ev$R2_ood, ev$R2_traj, el))
    }
  }
}
fwrite(rbindlist(job3_rows, fill = TRUE), "results/ablation_native/job3_nonunit_scale_expressivity.csv")
message("[OK] Saved results/ablation_native/job3_nonunit_scale_expressivity.csv")
}

# ==============================================================================
# JOB 4: Zero-Leakage Authentic Empirical Suite
#   1. auth_lynx_hare_1900_1920 (raw 21 yrs split FIRST, then spline strictly within train/OOD)
#   2. auth_huffaker_1963 (gauseR::huffaker_1963 60w real lab predator-prey mites)
#   3. auth_gause_1934_mixture (gauseR::gause_1934_book_f22 real Paramecium mixture)
#   4. real03_schmidt_pendulum & real05_yeast_glycolysis (zero pre-split leakage)
# ==============================================================================
if (has_job("4")) {
message("=== [JOB 4] Zero-Leakage Authentic Real-World Empirical Benchmarks ===")

build_leakfree_auth_lynx_hare <- function(fold_id = 1L, seed = 1900L, diff_method = "sgolay_spline") {
  df_raw <- read.csv("data/real_world/lynx_hare_annual_raw_1900_1920.csv", stringsAsFactors = FALSE)
  t_raw <- as.numeric(df_raw$time) # 0..20
  X_raw <- as.matrix(df_raw[, c("x1", "x2")])
  n_raw <- nrow(df_raw) # 21

  set.seed(seed)
  frac <- 0.65 + 0.10 * ((fold_id - 1L) %% 5L) / 4.0
  k_split <- max(12L, min(n_raw - 5L, floor(frac * n_raw)))

  # Split RAW years FIRST, then spline-upsample strictly within train vs OOD
  t_tr_raw <- t_raw[1:k_split]; X_tr_raw <- X_raw[1:k_split, , drop = FALSE]
  t_ood_raw <- t_raw[(k_split + 1L):n_raw]; X_ood_raw <- X_raw[(k_split + 1L):n_raw, , drop = FALSE]

  t_tr_grid <- seq(t_tr_raw[1], t_tr_raw[length(t_tr_raw)], by = 0.25)
  t_ood_grid <- seq(t_ood_raw[1], t_ood_raw[length(t_ood_raw)], by = 0.25)

  X_tr_dense <- cbind(
    x1 = pmax(1.0, stats::spline(t_tr_raw, X_tr_raw[, 1], xout = t_tr_grid, method = "natural")$y),
    x2 = pmax(1.0, stats::spline(t_tr_raw, X_tr_raw[, 2], xout = t_tr_grid, method = "natural")$y)
  )
  X_ood_dense <- cbind(
    x1 = pmax(1.0, stats::spline(t_ood_raw, X_ood_raw[, 1], xout = t_ood_grid, method = "natural")$y),
    x2 = pmax(1.0, stats::spline(t_ood_raw, X_ood_raw[, 2], xout = t_ood_grid, method = "natural")$y)
  )

  diff_fn <- if (diff_method == "sgolay_spline") estimate_derivatives_sgolay_spline else estimate_derivatives_tvregdiff
  Y_tr_dense <- diff_fn(X_tr_dense, t_tr_grid)
  Y_ood_dense <- diff_fn(X_ood_dense, t_ood_grid)

  n_tr <- nrow(X_tr_dense)
  idx_tr <- sort(sample(seq_len(n_tr), size = floor(0.85 * n_tr), replace = FALSE))
  idx_in <- setdiff(seq_len(n_tr), idx_tr)
  if (length(idx_in) < 3L) idx_in <- idx_tr

  list(
    system = list(
      sys_id = 101L, id = "auth_lynx_hare_1900_1920", name = "Authentic 1900-1920 Lynx-Hare (Zero-Leakage)",
      dim = 2L, state_names = c("x1", "x2"), feature_names = c("x1", "x2"), input_names = character(0),
      ground_truth = list(c("x1", "x1*x2"), c("x2", "x1*x2")), oracle_extra_basis = character(0)
    ),
    seed = seed, N_train = length(idx_tr),
    t_train = t_tr_grid[idx_tr], t_seq = t_tr_grid, X_seq = X_tr_dense,
    X_train = X_tr_dense[idx_tr, , drop = FALSE], Y_train = Y_tr_dense[idx_tr, , drop = FALSE],
    Y_train_clean = Y_tr_dense[idx_tr, , drop = FALSE],
    X_test_in = X_tr_dense[idx_in, , drop = FALSE], Y_test_in_clean = Y_tr_dense[idx_in, , drop = FALSE],
    Y_test_in_noisy = Y_tr_dense[idx_in, , drop = FALSE],
    X_test_ood = X_ood_dense, Y_test_ood_clean = Y_ood_dense, Y_test_ood_noisy = Y_ood_dense,
    rollout = list(times = t_ood_grid - t_ood_grid[1], X_true = X_ood_dense, x0 = X_ood_dense[1, ])
  )
}

build_leakfree_csv_task <- function(csv_path, ds_id, ds_name, gt_list, oracle_extra = character(0), fold_id = 1L, seed = 2000L) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  times <- as.numeric(df$time)
  X_states <- as.matrix(df[, c("x1", "x2")])
  n_tot <- nrow(X_states)
  set.seed(seed)
  frac <- 0.65 + 0.10 * ((fold_id - 1L) %% 5L) / 4.0
  n_tr_win <- max(14L, min(n_tot - 6L, floor(frac * n_tot)))
  idx_pool <- seq_len(n_tr_win)
  idx_ood <- seq(n_tr_win + 1L, n_tot)

  # Differentiate train & OOD strictly separately
  Y_tr_win <- estimate_derivatives_sgolay_spline(X_states[idx_pool, , drop = FALSE], times[idx_pool])
  Y_ood_win <- estimate_derivatives_sgolay_spline(X_states[idx_ood, , drop = FALSE], times[idx_ood])

  n_sub <- max(12L, floor(0.85 * length(idx_pool)))
  idx_tr_loc <- sort(sample(seq_along(idx_pool), size = n_sub, replace = FALSE))
  idx_in_loc <- setdiff(seq_along(idx_pool), idx_tr_loc)
  if (length(idx_in_loc) < 3L) idx_in_loc <- idx_tr_loc

  list(
    system = list(
      sys_id = 200L, id = ds_id, name = ds_name,
      dim = 2L, state_names = c("x1", "x2"), feature_names = c("x1", "x2"), input_names = character(0),
      ground_truth = gt_list, oracle_extra_basis = oracle_extra
    ),
    seed = seed, N_train = length(idx_tr_loc),
    t_train = times[idx_pool[idx_tr_loc]], t_seq = times[idx_pool], X_seq = X_states[idx_pool, , drop = FALSE],
    X_train = X_states[idx_pool[idx_tr_loc], , drop = FALSE], Y_train = Y_tr_win[idx_tr_loc, , drop = FALSE],
    Y_train_clean = Y_tr_win[idx_tr_loc, , drop = FALSE],
    X_test_in = X_states[idx_pool[idx_in_loc], , drop = FALSE],
    Y_test_in_clean = Y_tr_win[idx_in_loc, , drop = FALSE], Y_test_in_noisy = Y_tr_win[idx_in_loc, , drop = FALSE],
    X_test_ood = X_states[idx_ood, , drop = FALSE], Y_test_ood_clean = Y_ood_win, Y_test_ood_noisy = Y_ood_win,
    rollout = list(times = times[idx_ood] - times[idx_ood[1]], X_true = X_states[idx_ood, , drop = FALSE], x0 = X_states[idx_ood[1], ])
  )
}

j4_folds <- if (smoke_mode) 1L else 1:3
j4_methods <- c(
  "BGNLM_Full_Tuned", "SINDy_STLSQ_Poly", "UQ_SINDy_Poly",
  "UQ_SINDy_Bayes_Poly", "WSINDy_Poly", "SR_PySR", "SR_Operon"
)
job4_rows <- list()

dispatch_j4 <- function(mid, td) {
  if (mid == "BGNLM_Full_Tuned") {
    run_fbms_method(mid, td, cores = cores_arg, fast_mode = smoke_mode)
  } else if (mid %in% c("SINDy_STLSQ_Poly", "UQ_SINDy_Poly")) {
    run_sindy_baseline(mid, td, n_boot = if (smoke_mode) 25L else 80L)
  } else if (mid == "UQ_SINDy_Bayes_Poly") {
    run_uq_sindy_bayes(td, dict_type = "poly")
  } else if (mid == "WSINDy_Poly") {
    run_wsindy_baseline(td, dict_type = "poly")
  } else if (mid %in% c("SR_PySR", "SR_Operon")) {
    run_sr_baseline(mid, td, strict_native = TRUE)
  } else stop("Unknown method: ", mid)
}

j4_datasets <- c("auth_lynx_hare", "auth_huffaker_1963", "auth_gause_mixture", "real03_schmidt_pendulum", "real05_yeast_glycolysis")

for (ds in j4_datasets) {
  for (fid in j4_folds) {
    td <- if (ds == "auth_lynx_hare") {
      build_leakfree_auth_lynx_hare(fold_id = fid, seed = 1900L + fid)
    } else if (ds == "auth_huffaker_1963") {
      build_leakfree_csv_task("data/real_world/huffaker_predator_prey_authentic.csv",
        "auth_huffaker_1963", "Authentic Huffaker 1963 Mite Predator-Prey",
        list(c("x1", "x1*x2"), c("x2", "x1*x2")), fold_id = fid, seed = 1963L + fid)
    } else if (ds == "auth_gause_mixture") {
      build_leakfree_csv_task("data/real_world/gause_paramecium_mixture_authentic.csv",
        "auth_gause_mixture", "Authentic Gause 1934 Paramecium Competition",
        list(c("x1", "x1^2", "x1*x2"), c("x2", "x2^2", "x1*x2")), fold_id = fid, seed = 1934L + fid)
    } else {
      generate_real_task_data(ds, diff_method = "sgolay_spline", fold_id = fid, seed = 5000L + fid)
    }

    for (mid in j4_methods) {
      t0 <- proc.time()[3]
      fr <- run_safe_call(function() dispatch_j4(mid, td))
      el <- as.numeric(proc.time()[3] - t0)
      ev <- score_fit_full(fr, td, mid)
      job4_rows[[length(job4_rows) + 1L]] <- data.table(
        Dataset = ds, Fold_ID = fid, Method = mid, Engine = ev$Engine, Status = ev$Status,
        F1 = ev$F1, TPR = ev$TPR, FDR = ev$FDR, OC = ev$OC,
        R2_ood = ev$R2_ood, R2_traj = ev$R2_traj, Coverage_95_ood = ev$Coverage_95_ood,
        Runtime_s = el, Discovered_Atoms = ev$Discovered_Atoms
      )
      message(sprintf("[JOB 4] %s fold=%d | %s (%s) -> F1=%.3f FDR=%.3f R2_traj=%.3f (%.1fs)",
                      ds, fid, mid, ev$Engine, ev$F1, ifelse(is.na(ev$FDR), -1, ev$FDR), ev$R2_traj, el))
    }
  }
}
fwrite(rbindlist(job4_rows, fill = TRUE), "results/ablation_native/job4_authentic_real_world.csv")
message("[OK] Saved results/ablation_native/job4_authentic_real_world.csv")
}
message("=== ALL REQUESTED ABLATION JOBS COMPLETE ===")
