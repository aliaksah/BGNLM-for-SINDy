# ==============================================================================
# R/09_plotting_tables.R
# ICML Publication-Ready ggplot2 Figures & LaTeX Booktabs Tables
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

#' Format p-value for LaTeX scientific tables
format_pval_latex <- function(p) {
  ifelse(is.na(p), "--",
    ifelse(p < 0.0001, "< 10^{-4}",
      sprintf("%.4f", p)))
}

#' Clean term label for publication LaTeX tables
clean_term_for_latex <- function(term_str) {
  s <- gsub("^Method", "", as.character(term_str))
  s <- gsub("SystemGroupGroup_", "Group_", s, fixed = TRUE)
  s <- gsub("DiffMethodtvregdiff", "Diff: tvregdiff", s, fixed = TRUE)
  if (exists("relabel_method_str", mode = "function")) {
    s <- relabel_method_str(s)
  }
  gsub("_", "\\\\_", s)
}

#' Export Fixed-Effects Mixed-Model Contrasts to LaTeX Booktabs Table
export_mixed_effects_latex_table <- function(fixed_dt, out_path, caption_str, label_str) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(fixed_dt))
  if (!"df" %in% names(dt)) dt[, df := NA_real_]
  if (!"p.value" %in% names(dt)) dt[, p.value := NA_real_]
  if (!"p.dunnett" %in% names(dt)) dt[, p.dunnett := NA_real_]
  if (!"conf.low" %in% names(dt)) dt[, conf.low := estimate - 1.96 * std.error]
  if (!"conf.high" %in% names(dt)) dt[, conf.high := estimate + 1.96 * std.error]

  has_dunnett <- any(!is.na(dt$p.dunnett))
  col_spec <- if (has_dunnett) "llccccccc" else "llcccccc"
  hdr_extra <- if (has_dunnett) " & $p_{\\text{adj}}$" else ""

  lines <- c(
    "\\begin{table*}[t]",
    "\\centering",
    paste0("\\caption{", caption_str, "}"),
    paste0("\\label{", label_str, "}"),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.0pt}",
    "\\renewcommand{\\arraystretch}{0.86}",
    "\\resizebox{\\textwidth}{!}{%",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("Metric & Term / Contrast (vs.\\ \\texttt{BGNLM\\_Full\\_Tuned}) & $\\hat{\\beta}$ & SE & $df$ & $t$ & $95\\%$ CI & $p$-value", hdr_extra, " \\\\"),
    "\\midrule"
  )

  for (i in seq_len(nrow(dt))) {
    r <- dt[i]
    term_clean <- clean_term_for_latex(r$term)
    metric_clean <- gsub("_", "\\\\_", as.character(r$Metric))
    df_str <- if (is.na(r$df)) "--" else sprintf("%.1f", r$df)
    stat_str <- if (is.na(r$statistic)) "--" else sprintf("%.2f", r$statistic)
    ci_str <- sprintf("[%.3f, %.3f]", r$conf.low, r$conf.high)
    p_str <- paste0("$", format_pval_latex(r$p.value), "$")
    padj_col <- if (has_dunnett) {
      if (is.na(r$p.dunnett)) " & --" else paste0(" & $", format_pval_latex(r$p.dunnett), "$")
    } else ""

    lines <- c(lines, sprintf(
      "%s & \\texttt{%s} & %.3f & %.3f & %s & %s & %s & %s%s \\\\",
      metric_clean, term_clean, r$estimate, r$std.error, df_str, stat_str, ci_str, p_str, padj_col
    ))
  }

  lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table*}")
  writeLines(lines, out_path)
  invisible(out_path)
}

#' Export ICC Variance Decomposition Table to LaTeX
export_icc_latex_table <- function(icc_synth, icc_real, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  dt <- rbindlist(list(icc_synth, icc_real), fill = TRUE)
  has_sing <- "Singular" %in% names(dt)
  col_spec <- if (has_sing) "llcccc" else "llccc"
  sing_hdr <- if (has_sing) " & Singular" else ""
  lines <- c(
    "\\begin{table}[t]",
    "\\centering",
    "\\caption{Complete Random-Effects Variance Decomposition and Intraclass Correlation Coefficients (ICC).}",
    "\\label{tab:mixed_effects_icc}",
    "\\scriptsize",
    "\\resizebox{\\columnwidth}{!}{%",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("Model & Random Effect Group & $\\hat{\\sigma}^2$ & $\\hat{\\sigma}$ & ICC", sing_hdr, " \\\\"),
    "\\midrule"
  )
  for (i in seq_len(nrow(dt))) {
    r <- dt[i]
    m_clean <- gsub("_", "\\\\_", r$Model)
    g_clean <- gsub("_", "\\\\_", r$Group)
    s_col <- if (has_sing) paste0(" & ", r$Singular) else ""
    lines <- c(lines, sprintf(
      "%s & \\texttt{%s} & %.4f & %.4f & %.3f%s \\\\",
      m_clean, g_clean, r$Variance, r$StdDev, r$ICC, s_col
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}%", "}", "\\end{table}")
  writeLines(lines, out_path)
  invisible(out_path)
}

#' Generate Publication-Quality Mixed-Effects Forest Plot (Estimated Contrasts vs BGNLM_Full_Tuned)
plot_mixed_effects_forest <- function(fixed_dt, out_pdf, title_str = "Mixed-Effects Contrasts Relative to BGNLM_Full_Tuned") {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(fixed_dt))
  dt <- dt[grepl("^Method", term)]
  if (nrow(dt) == 0L) return(invisible(NULL))

  dt[, Method_Label := gsub("^Method", "", term)]
  if (exists("relabel_method_str", mode = "function")) {
    dt[, Method_Label := relabel_method_str(Method_Label)]
  }
  if (!"conf.low" %in% names(dt)) dt[, conf.low := estimate - 1.96 * std.error]
  if (!"conf.high" %in% names(dt)) dt[, conf.high := estimate + 1.96 * std.error]

  p <- ggplot(dt, aes(x = estimate, y = reorder(Method_Label, estimate), color = Metric)) +
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "black", linewidth = 0.7) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), orientation = "y", width = 0.25, linewidth = 0.7) +
    geom_point(size = 2.4) +
    facet_wrap(~ Metric, scales = "free_x") +
    labs(
      title = title_str,
      subtitle = "Dashed zero line = BGNLM_Full_Tuned reference (95% CI with System/Dataset x Method random interaction)",
      x = "Fixed-Effect Contrast (Method - BGNLM_Full_Tuned)",
      y = "Competing Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "#f0f4f8"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  ggsave(out_pdf, p, width = 10, height = 6, bg = "white")
  ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 10, height = 6, dpi = 200, bg = "white")
  invisible(p)
}

#' Generate Summary Performance Figure across Methods & System Groups
plot_benchmark_summary <- function(df_synth, out_pdf = "results/figures/fig_synth_benchmark_overview.pdf") {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(df_synth))
  if (nrow(dt) == 0L) return(invisible(NULL))
  if (exists("relabel_method_str", mode = "function")) {
    dt[, Method := relabel_method_str(Method)]
  }

  agg <- dt[, .(
    Mean_F1 = mean(F1, na.rm = TRUE),
    SE_F1 = sd(F1, na.rm = TRUE) / sqrt(max(1L, .N)),
    Mean_FDR = mean(FDR, na.rm = TRUE),
    Mean_R2_ood = mean(R2_test_ood, na.rm = TRUE)
  ), by = .(Method, SystemGroup)]

  p <- ggplot(agg, aes(x = reorder(Method, Mean_F1), y = Mean_F1, fill = SystemGroup)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(
      aes(ymin = pmax(0, Mean_F1 - 1.96 * SE_F1), ymax = pmin(1, Mean_F1 + 1.96 * SE_F1)),
      position = position_dodge(width = 0.8), width = 0.25
    ) +
    coord_flip() +
    labs(
      title = "Structural Identification F1 Score Across Benchmark Regimes",
      x = "Method",
      y = "Mean F1 Score (+/- 1.96 SE)",
      fill = "System Regime"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top")

  ggsave(out_pdf, p, width = 9, height = 5.5, bg = "white")
  ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 9, height = 5.5, dpi = 200, bg = "white")
  invisible(p)
}

