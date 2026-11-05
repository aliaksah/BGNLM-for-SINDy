#!/usr/bin/env Rscript
# ==============================================================================
# scripts/run_mixed_models.R
# Runs LMM/GLMM Statistical Evaluation & Exports LaTeX Tables and Forest Plots
# ==============================================================================

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

  export_mixed_effects_latex_table(
    res_synth$fixed_effects[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")],
    out_path = "results/tables/table_mixed_effects_synth.tex",
    caption_str = "Linear Mixed-Effects Model (LMM) Contrasts on Synthetic Benchmark Suite (Reference: \\texttt{BGNLM\\_Full\\_Tuned}).",
    label_str = "tab:mixed_effects_synth"
  )

  if (nrow(res_synth$regime_effects) > 0L) {
    export_mixed_effects_latex_table(
      res_synth$regime_effects,
      out_path = "results/tables/table_mixed_effects_regime.tex",
      caption_str = "Regime Moderation Interaction LMM (\\texttt{Method * SystemGroup}) Testing Polynomial vs. Non-Polynomial Dynamics.",
      label_str = "tab:mixed_effects_regime"
    )
  }

  plot_mixed_effects_forest(
    res_synth$fixed_effects[Metric %in% c("F1", "FDR", "R2_test_ood", "R2_traj")],
    out_pdf = "results/figures/fig_mixed_effects_forest_synth.pdf",
    title_str = "Experiment 1 (Synthetic ODEs): Mixed-Effects Contrasts vs. BGNLM_Full_Tuned"
  )
  plot_benchmark_summary(df_synth, out_pdf = "results/figures/fig_synth_benchmark_overview.pdf")
  message("[OK] Exported synthetic mixed-effects tables and forest plots.")
}

if (file.exists(real_file)) {
  df_real <- readRDS(real_file)
  message("[MIXED MODELS] Analyzing real-world empirical dataset results (N = ", nrow(df_real), " rows)...")
  res_real <- analyze_real_mixed_effects(df_real)
  icc_real <- res_real$icc_table

  fwrite(res_real$fixed_effects, "results/aggregated/mixed_effects_real_fixed.csv")
  if (nrow(res_real$dunnett_contrasts) > 0L) {
    fwrite(res_real$dunnett_contrasts, "results/aggregated/mixed_effects_real_dunnett.csv")
  }

  export_mixed_effects_latex_table(
    res_real$fixed_effects,
    out_path = "results/tables/table_mixed_effects_real.tex",
    caption_str = "Linear Mixed-Effects Model (LMM) Contrasts on Real-World Empirical Datasets (Reference: \\texttt{BGNLM\\_Full\\_Tuned}).",
    label_str = "tab:mixed_effects_real"
  )

  plot_mixed_effects_forest(
    res_real$fixed_effects,
    out_pdf = "results/figures/fig_mixed_effects_forest_real.pdf",
    title_str = "Experiment 2 (Real-World Datasets): Mixed-Effects Contrasts vs. BGNLM_Full_Tuned"
  )
  message("[OK] Exported real-data mixed-effects tables and forest plots.")
}

if (nrow(icc_synth) > 0L || nrow(icc_real) > 0L) {
  export_icc_latex_table(icc_synth, icc_real, "results/tables/table_mixed_effects_icc.tex")
}

