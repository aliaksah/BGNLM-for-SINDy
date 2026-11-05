# ==============================================================================
# R/08_mixed_effects_eval.R
# Mixed-Effects Statistical Evaluation Framework (lme4, lmerTest, glmmTMB, emmeans)
# ==============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(emmeans)
  library(broom.mixed)
  library(data.table)
})

#' Ensure Method factor has BGNLM_Full_Tuned as reference baseline level
prepare_mixed_model_data <- function(df, ref_method = "BGNLM_Full_Tuned") {
  dt <- as.data.table(copy(df))
  all_methods <- unique(as.character(dt$Method))
  if (ref_method %in% all_methods) {
    dt[, Method := relevel(factor(Method), ref = ref_method)]
  } else {
    dt[, Method := factor(Method)]
  }
  dt
}

#' Compute Intraclass Correlation Coefficients (ICC) from a fitted lme4 model
extract_random_effects_icc <- function(fit_mer, model_label = "LMM") {
  vc <- as.data.frame(lme4::VarCorr(fit_mer))
  total_var <- sum(vc$vcov, na.rm = TRUE)
  if (!is.finite(total_var) || total_var <= 0) total_var <- 1.0
  data.table(
    Model = model_label,
    Group = vc$grp,
    Variance = vc$vcov,
    StdDev = vc$sdcor,
    ICC = vc$vcov / total_var
  )
}

#' Fit Linear Mixed-Effects Model (LMM) cleanly with automatic random-effects fallback for small smoke-test grids
safe_fit_lmm <- function(formula_full, formula_simple, data) {
  fit <- tryCatch(
    suppressWarnings(lmerTest::lmer(formula_full, data = data, REML = TRUE)),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      suppressWarnings(lmerTest::lmer(formula_simple, data = data, REML = TRUE)),
      error = function(e) lm(nobars(formula_simple), data = data)
    )
  }
  fit
}

#' Run Full Mixed-Effects Analysis for Experiment 1 (Synthetic Benchmark Suite)
analyze_synthetic_mixed_effects <- function(df_synth) {
  dt <- prepare_mixed_model_data(df_synth, ref_method = "BGNLM_Full_Tuned")
  dt[, Seed_ID := paste0(System, "_N", SampleSize, "_eta", Noise, "_s", Seed)]
  dt[, log10_Noise := log10(pmax(as.numeric(Noise), 1e-6))]
  dt[, log10_N := log10(pmax(as.numeric(SampleSize), 10))]
  dt[, One_Minus_FDR := 1.0 - FDR]

  # Check if log10_Noise / log10_N vary in the dataset (smoke test may have 1 noise & 1 N)
  covars <- character(0)
  if (uniqueN(dt$log10_Noise) > 1L) covars <- c(covars, "log10_Noise")
  if (uniqueN(dt$log10_N) > 1L) covars <- c(covars, "log10_N")
  covar_str <- if (length(covars) > 0L) paste("+", paste(covars, collapse = " + ")) else ""

  metrics <- c("F1", "TPR", "FDR", "R2_test_in", "R2_test_ood", "R2_traj")
  lmm_fits <- list()
  tidy_rows <- list()
  dunnett_rows <- list()
  icc_rows <- list()

  for (m_name in metrics) {
    if (!m_name %in% names(dt)) next
    f_full <- as.formula(paste0(m_name, " ~ Method ", covar_str, " + (1 | System) + (1 | System:Seed_ID)"))
    f_simp <- as.formula(paste0(m_name, " ~ Method ", covar_str, " + (1 | Seed_ID)"))
    fit_m <- safe_fit_lmm(f_full, f_simp, dt)
    lmm_fits[[m_name]] <- fit_m

    td <- as.data.table(broom.mixed::tidy(fit_m, effects = "fixed", conf.int = TRUE))
    td[, Metric := m_name]
    tidy_rows[[length(tidy_rows) + 1L]] <- td

    if (inherits(fit_m, "merMod")) {
      icc_rows[[length(icc_rows) + 1L]] <- extract_random_effects_icc(fit_m, model_label = paste0("Synth_", m_name))
    }

    # Dunnett-adjusted treatment vs. control (BGNLM_Full_Tuned) contrasts via emmeans
    emm <- tryCatch({
      em_obj <- suppressMessages(emmeans::emmeans(fit_m, specs = "Method"))
      contr <- emmeans::contrast(em_obj, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      as.data.table(summary(contr, infer = c(TRUE, TRUE)))
    }, error = function(e) NULL)
    if (!is.null(emm)) {
      emm[, Metric := m_name]
      dunnett_rows[[length(dunnett_rows) + 1L]] <- emm
    }
  }

  # Regime Moderation Interaction Model: F1 ~ Method * SystemGroup
  fit_regime <- NULL
  tidy_regime <- data.table()
  if ("SystemGroup" %in% names(dt) && uniqueN(dt$SystemGroup) > 1L) {
    dt[, SystemGroup := factor(SystemGroup)]
    f_reg_full <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | System) + (1 | System:Seed_ID)"))
    f_reg_simp <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | Seed_ID)"))
    fit_regime <- safe_fit_lmm(f_reg_full, f_reg_simp, dt)
    tidy_regime <- as.data.table(broom.mixed::tidy(fit_regime, effects = "fixed", conf.int = TRUE))
    tidy_regime[, Metric := "F1_Regime_Interaction"]
  }

  # Binomial GLMM for Exact Model Recovery (I_exact in {0, 1})
  fit_glmm <- NULL
  tidy_glmm <- data.table()
  if ("Exact_Recovery" %in% names(dt) && uniqueN(dt$Exact_Recovery) > 1L) {
    f_glmm <- as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | Seed_ID)"))
    fit_glmm <- tryCatch(
      suppressWarnings(glmmTMB::glmmTMB(f_glmm, data = dt, family = binomial(link = "logit"))),
      error = function(e) {
        tryCatch(
          suppressWarnings(lme4::glmer(f_glmm, data = dt, family = binomial(link = "logit"))),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(fit_glmm)) {
      tidy_glmm <- tryCatch(
        as.data.table(broom.mixed::tidy(fit_glmm, effects = "fixed", conf.int = TRUE)),
        error = function(e) data.table()
      )
      if (nrow(tidy_glmm) > 0L) tidy_glmm[, Metric := "Exact_Recovery_Logit"]
    }
  }

  list(
    lmm_fits = lmm_fits,
    fit_regime = fit_regime,
    fit_glmm = fit_glmm,
    fixed_effects = rbindlist(tidy_rows, fill = TRUE),
    regime_effects = tidy_regime,
    glmm_effects = tidy_glmm,
    dunnett_contrasts = rbindlist(dunnett_rows, fill = TRUE),
    icc_table = rbindlist(icc_rows, fill = TRUE)
  )
}

#' Run Full Mixed-Effects Analysis for Experiment 2 (Real-World Empirical Datasets)
analyze_real_mixed_effects <- function(df_real) {
  dt <- prepare_mixed_model_data(df_real, ref_method = "BGNLM_Full_Tuned")
  dt[, Fold_UID := paste0(Dataset, "_", DiffMethod, "_fold", Fold_ID)]

  has_multi_diff <- ("DiffMethod" %in% names(dt)) && (uniqueN(dt$DiffMethod) > 1L)
  diff_str <- if (has_multi_diff) "+ DiffMethod" else ""

  metrics <- c("R2_traj", "NRMSE_traj", "R2_test_ood", "F1")
  lmm_fits <- list()
  tidy_rows <- list()
  dunnett_rows <- list()
  icc_rows <- list()

  for (m_name in metrics) {
    if (!m_name %in% names(dt)) next
    f_full <- as.formula(paste0(m_name, " ~ Method ", diff_str, " + (1 | Dataset) + (1 | Dataset:Fold_UID)"))
    f_simp <- as.formula(paste0(m_name, " ~ Method ", diff_str, " + (1 | Fold_UID)"))
    fit_m <- safe_fit_lmm(f_full, f_simp, dt)
    lmm_fits[[m_name]] <- fit_m

    td <- as.data.table(broom.mixed::tidy(fit_m, effects = "fixed", conf.int = TRUE))
    td[, Metric := m_name]
    tidy_rows[[length(tidy_rows) + 1L]] <- td

    if (inherits(fit_m, "merMod")) {
      icc_rows[[length(icc_rows) + 1L]] <- extract_random_effects_icc(fit_m, model_label = paste0("Real_", m_name))
    }

    emm <- tryCatch({
      em_obj <- suppressMessages(emmeans::emmeans(fit_m, specs = "Method"))
      contr <- emmeans::contrast(em_obj, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      as.data.table(summary(contr, infer = c(TRUE, TRUE)))
    }, error = function(e) NULL)
    if (!is.null(emm)) {
      emm[, Metric := m_name]
      dunnett_rows[[length(dunnett_rows) + 1L]] <- emm
    }
  }

  list(
    lmm_fits = lmm_fits,
    fixed_effects = rbindlist(tidy_rows, fill = TRUE),
    dunnett_contrasts = rbindlist(dunnett_rows, fill = TRUE),
    icc_table = rbindlist(icc_rows, fill = TRUE)
  )
}

