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

#' Map raw checkpoint Method IDs to transparent publication labels
relabel_method_str <- function(x) {
  x <- as.character(x)
  map <- c(
    "SR_PySR" = "SR_Pool_Pareto",
    "SR_Operon" = "SR_Pool_MDL",
    "UQ_SINDy_Oracle" = "Bag_STLSQ_Oracle",
    "UQ_SINDy_Poly" = "Bag_STLSQ_Poly",
    "SINDy_PI_Rational" = "SINDy_RatPoly2"
  )
  for (k in names(map)) {
    x <- gsub(k, map[[k]], x, fixed = TRUE)
  }
  x
}

#' Compute Intraclass Correlation Coefficients (ICC) from a fitted lme4 model
extract_random_effects_icc <- function(fit_mer, model_label = "LMM") {
  vc <- as.data.frame(lme4::VarCorr(fit_mer))
  total_var <- sum(vc$vcov, na.rm = TRUE)
  if (!is.finite(total_var) || total_var <= 0) total_var <- 1.0
  is_sing <- lme4::isSingular(fit_mer)
  data.table(
    Model = model_label,
    Group = vc$grp,
    Variance = vc$vcov,
    StdDev = vc$sdcor,
    ICC = vc$vcov / total_var,
    Singular = ifelse(is_sing, "Yes", "No")
  )
}

#' Merge Dunnett-adjusted p-values onto tidy fixed-effects table
attach_dunnett_pvals <- function(td, emm) {
  td[, p.dunnett := NA_real_]
  if (is.null(emm) || nrow(emm) == 0L) return(td)
  for (i in seq_len(nrow(td))) {
    trm <- as.character(td$term[i])
    if (startsWith(trm, "Method")) {
      m_lvl <- sub("^Method", "", trm)
      hit <- emm[grepl(paste0("^", m_lvl, "\\b"), contrast) | contrast == paste0(m_lvl, " - BGNLM_Full_Tuned")]
      if (nrow(hit) >= 1L && "p.value" %in% names(hit)) {
        td[i, p.dunnett := hit$p.value[1L]]
      }
    }
  }
  td
}

#' Fit Linear Mixed-Effects Model (LMM) cleanly with random-effects fallback and singularity logging
safe_fit_lmm <- function(formula_full, formula_simple, data, label = "LMM") {
  fit <- tryCatch(
    lmerTest::lmer(formula_full, data = data, REML = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      lmerTest::lmer(formula_simple, data = data, REML = TRUE),
      error = function(e) lm(lme4::nobars(formula_simple), data = data)
    )
  }
  if (inherits(fit, "merMod") && lme4::isSingular(fit)) {
    message("[LMM NOTE] Singular fit detected in ", label, ": ", paste(deparse(formula(fit)), collapse = " "))
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
    f_full <- as.formula(paste0(m_name, " ~ Method ", covar_str, " + (1 | System) + (1 | System:Method) + (1 | System:Seed_ID)"))
    f_simp <- as.formula(paste0(m_name, " ~ Method ", covar_str, " + (1 | Seed_ID)"))
    fit_m <- safe_fit_lmm(f_full, f_simp, dt, label = paste0("Synth_", m_name))
    lmm_fits[[m_name]] <- fit_m

    td <- as.data.table(broom.mixed::tidy(fit_m, effects = "fixed", conf.int = TRUE))
    td[, Metric := m_name]

    if (inherits(fit_m, "merMod")) {
      icc_rows[[length(icc_rows) + 1L]] <- extract_random_effects_icc(fit_m, model_label = paste0("Synth_", m_name))
    }

    # Dunnett-adjusted treatment vs. control (BGNLM_Full_Tuned) contrasts via emmeans
    emm <- tryCatch({
      em_obj <- suppressMessages(emmeans::emmeans(fit_m, specs = "Method", lmer.df = "satterthwaite", lmerTest.limit = 150000))
      contr <- emmeans::contrast(em_obj, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      as.data.table(summary(contr, infer = c(TRUE, TRUE)))
    }, error = function(e) NULL)
    td <- attach_dunnett_pvals(td, emm)
    tidy_rows[[length(tidy_rows) + 1L]] <- td
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
    f_reg_full <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | System) + (1 | System:Method) + (1 | System:Seed_ID)"))
    f_reg_simp <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | Seed_ID)"))
    fit_regime <- safe_fit_lmm(f_reg_full, f_reg_simp, dt, label = "Synth_Regime")
    tidy_regime <- as.data.table(broom.mixed::tidy(fit_regime, effects = "fixed", conf.int = TRUE))
    tidy_regime[, Metric := "F1_Regime"]
    tidy_regime[, p.dunnett := p.adjust(p.value, method = "holm")]
  }

  # Binomial GLMM for Exact Model Recovery (I_exact in {0, 1})
  fit_glmm <- NULL
  tidy_glmm <- data.table()
  if ("Exact_Recovery" %in% names(dt) && uniqueN(dt$Exact_Recovery) > 1L) {
    f_glmm <- as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | System) + (1 | System:Method)"))
    fit_glmm <- tryCatch(
      glmmTMB::glmmTMB(f_glmm, data = dt, family = binomial(link = "logit")),
      error = function(e) {
        tryCatch(
          glmmTMB::glmmTMB(as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | Seed_ID)")), data = dt, family = binomial(link = "logit")),
          error = function(e2) NULL
        )
      }
    )
    if (!is.null(fit_glmm)) {
      tidy_glmm <- tryCatch(
        as.data.table(broom.mixed::tidy(fit_glmm, effects = "fixed", conf.int = TRUE)),
        error = function(e) data.table()
      )
      if (nrow(tidy_glmm) > 0L) {
        tidy_glmm[, Metric := "Exact_Recovery_Logit"]
        tidy_glmm[, p.dunnett := ifelse(startsWith(term, "Method"), p.adjust(p.value, method = "holm"), NA_real_)]
      }
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
    f_full <- as.formula(paste0(m_name, " ~ Method ", diff_str, " + (1 | Dataset) + (1 | Dataset:Method) + (1 | Dataset:Fold_UID)"))
    f_simp <- as.formula(paste0(m_name, " ~ Method ", diff_str, " + (1 | Fold_UID)"))
    fit_m <- safe_fit_lmm(f_full, f_simp, dt, label = paste0("Real_", m_name))
    lmm_fits[[m_name]] <- fit_m

    td <- as.data.table(broom.mixed::tidy(fit_m, effects = "fixed", conf.int = TRUE))
    td[, Metric := m_name]

    if (inherits(fit_m, "merMod")) {
      icc_rows[[length(icc_rows) + 1L]] <- extract_random_effects_icc(fit_m, model_label = paste0("Real_", m_name))
    }

    emm <- tryCatch({
      em_obj <- suppressMessages(emmeans::emmeans(fit_m, specs = "Method", lmer.df = "satterthwaite", lmerTest.limit = 10000))
      contr <- emmeans::contrast(em_obj, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      as.data.table(summary(contr, infer = c(TRUE, TRUE)))
    }, error = function(e) NULL)
    td <- attach_dunnett_pvals(td, emm)
    tidy_rows[[length(tidy_rows) + 1L]] <- td
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

