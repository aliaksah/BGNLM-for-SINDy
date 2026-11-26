#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_mixed_models.R
# Runs LMM/GLMM Statistical Evaluation & Exports LaTeX Tables and Forest Plots
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

source("R/08_mixed_effects_eval.R")
source("R/09_plotting_tables.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

synth_file <- "results/aggregated/synth_results.rds"
real_file <- "results/aggregated/real_results.rds"

icc_synth <- data.table()
icc_real <- data.table()

if (file.exists(synth_file)) {
  df_synth <- readRDS(synth_file)
  message("[MIXED MODELS] Analyzing synthetic benchmark results (N = ", nrow(df_synth), " rows)...")
  res_synth <- analyze_synthetic_mixed_effects(df_synth)
  icc_synth <- res_synth$icc_table

  fwrite(res_synth$fixed_effects, "results/aggregated/mixed_effects_synth_fixed.csv")
  if (nrow(res_synth$dunnett_contrasts) > 0L) {
    fwrite(res_synth$dunnett_contrasts, "results/aggregated/mixed_effects_synth_dunnett.csv")
  }
  if (nrow(res_synth$regime_effects) > 0L) {
    fwrite(res_synth$regime_effects, "results/aggregated/mixed_effects_synth_regime.csv")
  }
  if (!is.null(res_synth$regime_dunnett) && nrow(res_synth$regime_dunnett) > 0L) {
    fwrite(res_synth$regime_dunnett, "results/aggregated/mixed_effects_synth_regime_dunnett.csv")
    rd <- copy(res_synth$regime_dunnett)
    rd[, Metric := as.character(SystemGroup)]
    rd[, term := paste0("Method", sub(" - BGNLM_Full_Tuned$", "", contrast))]
    if ("SE" %in% names(rd)) setnames(rd, "SE", "std.error")
    rd[, statistic := t.ratio]
    rd[, p.dunnett := p.value]
    if (all(c("lower.CL", "upper.CL") %in% names(rd))) {
      setnames(rd, c("lower.CL", "upper.CL"), c("conf.low", "conf.high"))
    }
    export_mixed_effects_latex_table(
      rd,
      out_path = "results/tables/tbl_within_regime_contrasts.tex",
      caption_str = "Within-Regime Dunnett-Adjusted LMM Contrasts for $F_1$ (\\texttt{Method} vs.\\ \\texttt{BGNLM\\_Full\\_Tuned}) Conditioned on Dynamical Regime Group.",
      label_str = "tab:within_regime_contrasts"
    )
  }
  if (!is.null(res_synth$ordbeta_contrasts) && nrow(res_synth$ordbeta_contrasts) > 0L) {
    fwrite(res_synth$ordbeta_contrasts, "results/aggregated/mixed_effects_synth_ordbeta.csv")
    ob <- copy(res_synth$ordbeta_contrasts)
    ob[, term := paste0("Method", sub(" - BGNLM_Full_Tuned$", "", contrast))]
    if ("SE" %in% names(ob)) setnames(ob, "SE", "std.error")
    ob[, statistic := if ("z.ratio" %in% names(ob)) z.ratio else t.ratio]
    ob[, p.dunnett := p.value]
    if (all(c("asymp.LCL", "asymp.UCL") %in% names(ob))) {
      setnames(ob, c("asymp.LCL", "asymp.UCL"), c("conf.low", "conf.high"))
    } else if (all(c("lower.CL", "upper.CL") %in% names(ob))) {
      setnames(ob, c("lower.CL", "upper.CL"), c("conf.low", "conf.high"))
    }
    export_mixed_effects_latex_table(
      ob,
      out_path = "results/tables/tbl_ordbeta_sensitivity.tex",
      caption_str = "Ordered Beta GLMM Sensitivity Check (\\texttt{glmmTMB::ordbeta}, Link-Scale Logit Contrasts vs.\\ \\texttt{BGNLM\\_Full\\_Tuned}) Confirming Full Concordance on $[0,1]$-Bounded Outcomes.",
      label_str = "tab:ordbeta_sensitivity"
    )
  }
  if (!is.null(res_synth$glmm_effects) && nrow(res_synth$glmm_effects) > 0L) {
    fwrite(res_synth$glmm_effects, "results/aggregated/mixed_effects_synth_glmm.csv")
    export_mixed_effects_latex_table(
      res_synth$glmm_effects,
      out_path = "results/tables/table_mixed_effects_glmm.tex",
      caption_str = "Binomial GLMM Log-Odds Contrasts for Exact System Recovery $I_{\\text{exact}} \\in \\{0,1\\}$ (Reference: \\texttt{BGNLM\\_Full\\_Tuned}).",
      label_str = "tab:mixed_effects_glmm"
    )
  }

  # Held-out 9-system validation LMM (excluding tuning systems sys02, sys04, sys08)
  df_held9 <- df_synth[!System %in% c("sys02_lorenz63", "sys04_lotka_volterra", "sys08_rossler_lorenz_hybrid")]
  if (nrow(df_held9) > 0L) {
    res_held9 <- analyze_synthetic_mixed_effects(df_held9)
    fwrite(res_held9$fixed_effects[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")], "results/aggregated/mixed_effects_synth_heldout9.csv")
  }

  export_mixed_effects_latex_table(
    res_synth$fixed_effects[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")],
    out_path = "results/tables/table_mixed_effects_synth.tex",
    caption_str = "Linear Mixed-Effects Model (LMM) Contrasts on Synthetic Benchmark Suite with \\texttt{(1|System:Method)} Random Interaction (Reference: \\texttt{BGNLM\\_Full\\_Tuned}).",
    label_str = "tab:mixed_effects_synth"
  )

  if (nrow(res_synth$regime_effects) > 0L) {
    export_mixed_effects_latex_table(
      res_synth$regime_effects,
      out_path = "results/tables/table_mixed_effects_regime.tex",
      caption_str = "Regime Moderation Interaction LMM (\\texttt{Method * SystemGroup}) Testing Polynomial vs.\\ Non-Polynomial vs.\\ 5D Dynamics (\\texttt{Group\\_C} comprises 5D Lorenz-96).",
      label_str = "tab:mixed_effects_regime"
    )
  }

  plot_mixed_effects_forest(
    res_synth$fixed_effects[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")],
    dunnett_dt = res_synth$dunnett_contrasts[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")],
    out_pdf = "results/figures/fig_mixed_effects_forest_synth.pdf"
  )
  plot_benchmark_summary(df_synth, out_pdf = "results/figures/fig_synth_benchmark_overview.pdf")
  message("[OK] Exported synthetic mixed-effects tables and forest plots.")
}

if (file.exists(real_file)) {
  df_real <- readRDS(real_file)
  message("[MIXED MODELS] Analyzing quasi-empirical ODE dataset results (N = ", nrow(df_real), " rows)...")
  res_real <- analyze_real_mixed_effects(df_real)
  icc_real <- res_real$icc_table

  fwrite(res_real$fixed_effects, "results/aggregated/mixed_effects_real_fixed.csv")
  if (nrow(res_real$dunnett_contrasts) > 0L) {
    fwrite(res_real$dunnett_contrasts, "results/aggregated/mixed_effects_real_dunnett.csv")
  }

  export_mixed_effects_latex_table(
    res_real$fixed_effects,
    out_path = "results/tables/table_mixed_effects_real.tex",
    caption_str = "Linear Mixed-Effects Model (LMM) Contrasts on Quasi-Empirical Literature ODE Benchmarks with \\texttt{(1|Dataset:Method)} Random Interaction (Reference: \\texttt{BGNLM\\_Full\\_Tuned}).",
    label_str = "tab:mixed_effects_real"
  )

  plot_mixed_effects_forest(
    res_real$fixed_effects,
    dunnett_dt = res_real$dunnett_contrasts,
    out_pdf = "results/figures/fig_mixed_effects_forest_real.pdf"
  )
  message("[OK] Exported quasi-empirical mixed-effects tables and forest plots.")
}

if (nrow(icc_synth) > 0L || nrow(icc_real) > 0L) {
  export_icc_latex_table(icc_synth, icc_real, "results/tables/table_mixed_effects_icc.tex")
}

