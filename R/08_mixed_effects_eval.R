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

    # Dunnett-adjusted treatment vs. control (BGNLM_Full_Tuned) contrasts + simultaneous vs unadjusted 95% CIs
    emm <- tryCatch({
      em_obj <- suppressMessages(emmeans::emmeans(fit_m, specs = "Method", lmer.df = "satterthwaite", lmerTest.limit = 150000))
      contr <- emmeans::contrast(em_obj, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      s_dun <- as.data.table(summary(contr, infer = c(TRUE, TRUE), adjust = "dunnett"))
      s_non <- as.data.table(summary(contr, infer = c(TRUE, FALSE), adjust = "none"))
      if (all(c("lower.CL", "upper.CL") %in% names(s_dun))) {
        setnames(s_dun, c("lower.CL", "upper.CL"), c("dunnett.low", "dunnett.high"))
      }
      if (all(c("lower.CL", "upper.CL") %in% names(s_non)) && nrow(s_non) == nrow(s_dun)) {
        s_dun[, unadj.low := s_non$lower.CL]
        s_dun[, unadj.high := s_non$upper.CL]
      }
      s_dun
    }, error = function(e) NULL)
    td <- attach_dunnett_pvals(td, emm)
    tidy_rows[[length(tidy_rows) + 1L]] <- td
    if (!is.null(emm)) {
      emm[, Metric := m_name]
      dunnett_rows[[length(dunnett_rows) + 1L]] <- emm
    }
  }

  # Regime Moderation Interaction Model + Within-Regime Dunnett Contrasts: F1 ~ Method * SystemGroup
  fit_regime <- NULL
  tidy_regime <- data.table()
  regime_dunnett <- data.table()
  if ("SystemGroup" %in% names(dt) && uniqueN(dt$SystemGroup) > 1L) {
    dt[, SystemGroup := factor(SystemGroup)]
    f_reg_full <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | System) + (1 | System:Method) + (1 | System:Seed_ID)"))
    f_reg_simp <- as.formula(paste0("F1 ~ Method * SystemGroup ", covar_str, " + (1 | Seed_ID)"))
    fit_regime <- safe_fit_lmm(f_reg_full, f_reg_simp, dt, label = "Synth_Regime")
    tidy_regime <- as.data.table(broom.mixed::tidy(fit_regime, effects = "fixed", conf.int = TRUE))
    tidy_regime[, Metric := "F1_Regime"]
    tidy_regime[, p.dunnett := p.adjust(p.value, method = "holm")]

    regime_dunnett <- tryCatch({
      em_reg <- suppressMessages(emmeans::emmeans(
        fit_regime, specs = ~ Method | SystemGroup,
        lmer.df = "satterthwaite", lmerTest.limit = 150000
      ))
      contr_reg <- emmeans::contrast(em_reg, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
      as.data.table(summary(contr_reg, infer = c(TRUE, TRUE), adjust = "dunnett"))
    }, error = function(e) data.table())
  }

  # Binomial GLMM for Exact Model Recovery (I_exact in {0, 1}) with block RE (1 | System:Seed_ID)
  fit_glmm <- NULL
  tidy_glmm <- data.table()
  if ("Exact_Recovery" %in% names(dt) && uniqueN(dt$Exact_Recovery) > 1L) {
    glmm_ctrl <- glmmTMB::glmmTMBControl(parallel = min(4L, max(1L, parallel::detectCores(logical = FALSE))))
    f_glmm_full <- as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | System) + (1 | System:Method) + (1 | System:Seed_ID)"))
    f_glmm_red <- as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | System) + (1 | System:Seed_ID)"))
    fit_glmm <- tryCatch(
      glmmTMB::glmmTMB(f_glmm_full, data = dt, family = binomial(link = "logit"), control = glmm_ctrl),
      error = function(e) {
        tryCatch(
          glmmTMB::glmmTMB(f_glmm_red, data = dt, family = binomial(link = "logit"), control = glmm_ctrl),
          error = function(e2) {
            tryCatch(
              glmmTMB::glmmTMB(as.formula(paste0("Exact_Recovery ~ Method ", covar_str, " + (1 | System) + (1 | System:Method)")), data = dt, family = binomial(link = "logit"), control = glmm_ctrl),
              error = function(e3) NULL
            )
          }
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

  # Ordered Beta GLMM sensitivity check (Kubinec 2023) on [0,1]-bounded F1, FDR, R2_traj
  ordbeta_rows <- list()
  lmm_dun_dt <- rbindlist(dunnett_rows, fill = TRUE)
  ob_ctrl <- glmmTMB::glmmTMBControl(parallel = min(4L, max(1L, parallel::detectCores(logical = FALSE))))
  for (ob_m in c("F1", "FDR", "R2_traj")) {
    if (!ob_m %in% names(dt)) next
    dt_ob <- dt[is.finite(get(ob_m))]
    dt_ob[, y_ob := pmin(1.0, pmax(0.0, as.numeric(get(ob_m))))]
    # Fit cross-system random-interaction Ordered Beta GLMM (132 System:Method levels keep Laplace exact & fast)
    f_ob <- as.formula(paste0("y_ob ~ Method ", covar_str, " + (1 | System) + (1 | System:Method)"))
    fit_ob <- tryCatch(
      glmmTMB::glmmTMB(f_ob, data = dt_ob, family = glmmTMB::ordbeta(), control = ob_ctrl),
      error = function(e) {
        tryCatch(glmmTMB::glmmTMB(as.formula(paste0("y_ob ~ Method ", covar_str, " + (1 | System)")), data = dt_ob, family = glmmTMB::ordbeta(), control = ob_ctrl), error = function(e2) NULL)
      }
    )
    if (!is.null(fit_ob)) {
      em_ob <- tryCatch({
        em_o <- suppressMessages(emmeans::emmeans(fit_ob, specs = "Method"))
        contr_o <- emmeans::contrast(em_o, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett")
        dt_o <- as.data.table(summary(contr_o, infer = c(TRUE, TRUE), adjust = "dunnett"))
        dt_o[, Metric := ob_m]
        dt_o
      }, error = function(e) NULL)
      if (!is.null(em_ob)) ordbeta_rows[[length(ordbeta_rows) + 1L]] <- em_ob
    }
  }
  ordbeta_dt <- rbindlist(ordbeta_rows, fill = TRUE)

  list(
    lmm_fits = lmm_fits,
    fit_regime = fit_regime,
    fit_glmm = fit_glmm,
    fixed_effects = rbindlist(tidy_rows, fill = TRUE),
    regime_effects = tidy_regime,
    regime_dunnett = regime_dunnett,
    glmm_effects = tidy_glmm,
    dunnett_contrasts = lmm_dun_dt,
    ordbeta_contrasts = ordbeta_dt,
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
      s_dun <- as.data.table(summary(contr, infer = c(TRUE, TRUE), adjust = "dunnett"))
      s_non <- as.data.table(summary(contr, infer = c(TRUE, FALSE), adjust = "none"))
      if (all(c("lower.CL", "upper.CL") %in% names(s_dun))) {
        setnames(s_dun, c("lower.CL", "upper.CL"), c("dunnett.low", "dunnett.high"))
      }
      if (all(c("lower.CL", "upper.CL") %in% names(s_non)) && nrow(s_non) == nrow(s_dun)) {
        s_dun[, unadj.low := s_non$lower.CL]
        s_dun[, unadj.high := s_non$upper.CL]
      }
      s_dun
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

