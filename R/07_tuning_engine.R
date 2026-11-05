# ==============================================================================
# R/07_tuning_engine.R
# Systematic Hyperparameter & Prior Tuning Engine for FBMS (BGNLM / BGLM)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

for (f in c(
  "R/01_ode_simulators.R",
  "R/02_basis_builders.R",
  "R/05_symbolic_canonicalize.R",
  "R/06_metrics_consistency.R",
  "R/03_fbms_wrappers.R"
)) {
  if (file.exists(f)) source(f)
  else if (file.exists(file.path("../..", f))) source(file.path("../..", f))
}

#' Candidate Hyperparameter Configurations Comparing Parameter Priors, Complexity Priors, & Operator Probs
get_bgnlm_tuning_candidates <- function() {
  list(
    list(
      config_id = "C01_GPrior_rN1_Algebraic",
      label = "Zellner g-prior (g=max(N,p^2), unknown sigma^2), r=N^-1, p_proj=0.0, p_inter=0.55",
      beta_prior_type = "g-prior",
      g_rule = "max_n_p2",
      var = "unknown",
      r_rule = "N^(-1.0)",
      p_inter = 0.55, p_mod = 0.35, p_proj = 0.00, p_new = 0.10,
      pop_max = 20L, max_depth_D = 2L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C02_HyperGN_rN1_Algebraic",
      label = "Mixture Hyper-g/n prior (unknown sigma^2), r=N^-1, p_proj=0.0, p_inter=0.55",
      beta_prior_type = "hyper-g-n",
      g_rule = "hyper-g-n",
      var = "unknown",
      r_rule = "N^(-1.0)",
      p_inter = 0.55, p_mod = 0.35, p_proj = 0.00, p_new = 0.10,
      pop_max = 20L, max_depth_D = 2L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C03_Robust_tCCH_rN1_Algebraic",
      label = "Robust tCCH prior (unknown sigma^2), r=N^-1, p_proj=0.0, p_inter=0.55",
      beta_prior_type = "robust",
      g_rule = "robust",
      var = "unknown",
      r_rule = "N^(-1.0)",
      p_inter = 0.55, p_mod = 0.35, p_proj = 0.00, p_new = 0.10,
      pop_max = 20L, max_depth_D = 2L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C04_JeffreysBIC_UnknownVar_rN1",
      label = "Jeffreys-BIC (unknown sigma^2), r=N^-1, p_proj=0.0, p_inter=0.55",
      beta_prior_type = "Jeffreys-BIC",
      g_rule = "bic",
      var = "unknown",
      r_rule = "N^(-1.0)",
      p_inter = 0.55, p_mod = 0.35, p_proj = 0.00, p_new = 0.10,
      pop_max = 20L, max_depth_D = 2L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C05_GPrior_rN15_StrongPenalty",
      label = "Zellner g-prior (unknown sigma^2), r=N^-1.5, p_proj=0.0, p_inter=0.55",
      beta_prior_type = "g-prior",
      g_rule = "max_n_p2",
      var = "unknown",
      r_rule = "N^(-1.5)",
      p_inter = 0.55, p_mod = 0.35, p_proj = 0.00, p_new = 0.10,
      pop_max = 25L, max_depth_D = 2L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C06_GPrior_rExp4_LowProj",
      label = "Zellner g-prior (unknown sigma^2), r=exp(-4), p_proj=0.05, p_inter=0.50",
      beta_prior_type = "g-prior",
      g_rule = "max_n_p2",
      var = "unknown",
      r_rule = "exp(-4)",
      p_inter = 0.50, p_mod = 0.35, p_proj = 0.05, p_new = 0.10,
      pop_max = 20L, max_depth_D = 3L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_extended", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C07_GPrior_HighProj_Abalation",
      label = "Zellner g-prior (unknown sigma^2), r=N^-0.5, p_proj=0.20, p_inter=0.40",
      beta_prior_type = "g-prior",
      g_rule = "max_n_p2",
      var = "unknown",
      r_rule = "N^(-0.5)",
      p_inter = 0.40, p_mod = 0.30, p_proj = 0.20, p_new = 0.10,
      pop_max = 15L, max_depth_D = 3L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_core", energy_prune_tol = 0.025
    ),
    list(
      config_id = "C08_IWSM_Baseline_FixedVar1",
      label = "Untuned IWSM Baseline (Jeffreys fixed sigma^2=1, sin_deg/cos_deg, pop.max=15)",
      beta_prior_type = "Jeffreys-BIC",
      g_rule = "fixed_var1",
      var = 1.0,
      r_rule = "N^(-1.0)",
      p_inter = 0.40, p_mod = 0.40, p_proj = 0.10, p_new = 0.10,
      pop_max = 15L, max_depth_D = 5L, P_generations = 12L, N_mcmc = 500L, runs = 6L,
      transforms_tier = "transforms_iwsm", energy_prune_tol = 1e-4
    )
  )
}

#' Run BGNLM Calibration Suite and Select Optimal Configuration
execute_bgnlm_tuning <- function(cores = 2L, fast_calibration = TRUE) {
  candidates <- get_bgnlm_tuning_candidates()
  cal_systems <- c("sys02_lorenz63", "sys08_rossler_lorenz_hybrid", "sys04_lotka_volterra")
  cal_noises <- c(0.01, 0.25)
  cal_N <- if (fast_calibration) c(250L) else c(250L, 1000L)
  cal_seeds <- if (fast_calibration) c(101L) else c(101L, 102L)

  rows <- list()
  for (cand in candidates) {
    message("[TUNING] Evaluating configuration: ", cand$config_id, " (", cand$label, ")")
    for (sys_id in cal_systems) {
      for (n_obs in cal_N) {
        for (eta in cal_noises) {
          for (sd_val in cal_seeds) {
            task_data <- generate_synthetic_task_data(
              system_id = sys_id,
              N_train = n_obs,
              noise_eta = eta,
              noise_regime = "derivative_snr",
              seed = sd_val,
              T_max = 15.0
            )
            N_num <- as.numeric(n_obs)
            p_cov <- ncol(task_data$X_train)
            beta_p <- if (cand$beta_prior_type == "g-prior") {
              list(type = "g-prior", g = max(N_num, p_cov^2))
            } else if (cand$beta_prior_type == "Jeffreys-BIC") {
              list(type = "Jeffreys-BIC", var = cand$var)
            } else {
              list(type = cand$beta_prior_type)
            }
            r_val <- switch(
              cand$r_rule,
              "N^(-0.5)" = N_num^(-0.5),
              "N^(-1.0)" = N_num^(-1.0),
              "N^(-1.5)" = N_num^(-1.5),
              "exp(-2)"  = exp(-2.0),
              "exp(-4)"  = exp(-4.0),
              1.0 / N_num
            )
            override <- cand
            override$beta_prior <- beta_p
            override$model_prior <- list(r = r_val)

            method_call <- if (cand$config_id == "C08_IWSM_Baseline_FixedVar1") "BGNLM_IWSM_Untuned" else "BGNLM_Full_Tuned"
            t0 <- proc.time()[3]
            fit_res <- run_fbms_method(
              method_id = method_call,
              task_data = task_data,
              cores = cores,
              tuning_override = override,
              fast_mode = fast_calibration
            )
            elapsed <- proc.time()[3] - t0

            st_m <- compute_structure_metrics(fit_res$discovered, task_data$system$ground_truth)
            pip_brier <- compute_pip_brier_score(fit_res$pip_list, task_data$system$ground_truth)
            pr_ood <- compute_prediction_metrics(fit_res$pred_ood, task_data$Y_test_ood_clean, fit_res$lower_ood, fit_res$upper_ood)

            rows[[length(rows) + 1L]] <- data.table(
              config_id = cand$config_id,
              label = cand$label,
              beta_prior_type = cand$beta_prior_type,
              r_rule = cand$r_rule,
              p_proj = cand$p_proj,
              p_inter = cand$p_inter,
              system_id = sys_id,
              N_train = n_obs,
              noise_eta = eta,
              seed = sd_val,
              TPR = st_m$TPR,
              FDR = st_m$FDR,
              F1 = st_m$F1,
              Exact_Recovery = st_m$Exact_Recovery,
              SHD = st_m$SHD,
              PIP_Brier = pip_brier,
              R2_ood = pr_ood$R2_clipped,
              Runtime_s = elapsed
            )
          }
        }
      }
    }
  }

  dt_all <- rbindlist(rows)
  summary_dt <- dt_all[, .(
    Mean_F1 = mean(F1, na.rm = TRUE),
    Mean_TPR = mean(TPR, na.rm = TRUE),
    Mean_FDR = mean(FDR, na.rm = TRUE),
    Exact_Rate = mean(Exact_Recovery, na.rm = TRUE),
    Mean_SHD = mean(SHD, na.rm = TRUE),
    Mean_PIP_Brier = mean(PIP_Brier, na.rm = TRUE),
    Mean_R2_ood = mean(R2_ood, na.rm = TRUE),
    Mean_Runtime_s = mean(Runtime_s, na.rm = TRUE)
  ), by = .(config_id, label, beta_prior_type, r_rule, p_proj, p_inter)]

  # Composite selection score: 0.50 * Mean_F1 + 0.20 * Exact_Rate + 0.20 * Mean_R2_ood - 0.10 * Mean_PIP_Brier
  summary_dt[, Composite_Score := 0.50 * Mean_F1 + 0.20 * Exact_Rate + 0.20 * Mean_R2_ood - 0.10 * fifelse(is.na(Mean_PIP_Brier), 0.25, Mean_PIP_Brier)]
  setorder(summary_dt, -Composite_Score)

  dir.create("results/aggregated", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
  fwrite(dt_all, "results/aggregated/tuning_runs_raw.csv")
  fwrite(summary_dt, "results/aggregated/tuning_summary.csv")

  best_id <- summary_dt$config_id[1]
  best_cand <- candidates[[which(vapply(candidates, function(x) x$config_id == best_id, logical(1)))[1]]]
  yaml::write_yaml(best_cand, "config/best_bgnlm_params.yaml")
  message("[TUNING] Locked optimal configuration: ", best_cand$config_id, " -> config/best_bgnlm_params.yaml")

  # Export LaTeX booktabs table of tuning ablation
  tex_lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Systematic BGNLM Hyperparameter and Prior Tuning on the Calibration Suite.}",
    "\\label{tab:bgnlm_tuning}",
    "\\small",
    "\\begin{tabular}{llcccccc}",
    "\\toprule",
    "Config ID & Parameter Prior & $r$ Penalty & $p_{\\text{proj}}$ & F1 $\\uparrow$ & FDR $\\downarrow$ & PIP Brier $\\downarrow$ & $R^2_{\\text{OOD}}$ $\\uparrow$ \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary_dt))) {
    r <- summary_dt[i]
    cfg_clean <- gsub("_", "\\\\_", r$config_id)
    prior_clean <- gsub("_", "\\\\_", r$beta_prior_type)
    tex_lines <- c(tex_lines, sprintf(
      "%s & %s & \\texttt{%s} & %.2f & %.3f & %.3f & %.3f & %.3f \\\\",
      cfg_clean, prior_clean, r$r_rule, r$p_proj,
      r$Mean_F1, r$Mean_FDR, r$Mean_PIP_Brier, r$Mean_R2_ood
    ))
  }
  tex_lines <- c(tex_lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex_lines, "results/tables/table_bgnlm_tuning.tex")

  list(raw = dt_all, summary = summary_dt, best = best_cand)
}

