# ==============================================================================
# R/09_plotting_tables.R
# ICML Publication-Ready ggplot2 Figures & LaTeX Booktabs Tables (Okabe-Ito)
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

#' Map raw method strings to human-readable academic labels
human_method_label <- function(m_str) {
  m_str <- as.character(m_str)
  map <- c(
    "BGNLM_Full_Tuned"    = "BGNLM (Tuned)",
    "BGLM_Oracle_Basis"   = "BGLM (Oracle)",
    "BGNLM_Hybrid_Poly"   = "BGNLM-Hybrid",
    "BGLM_Poly_Basis"     = "BGLM (Poly)",
    "BGLM_Overcomplete"   = "BGLM (Overcomplete)",
    "BGNLM_IWSM_Untuned"  = "BGNLM (Untuned IWSM)",
    "SINDy_STLSQ_Oracle"  = "STLSQ (Oracle)",
    "SINDy_STLSQ_Poly"    = "STLSQ (Poly)",
    "UQ_SINDy_Oracle"     = "Bag-STLSQ (Oracle)",
    "Bag_STLSQ_Oracle"    = "Bag-STLSQ (Oracle)",
    "UQ_SINDy_Poly"       = "Bag-STLSQ (Poly)",
    "Bag_STLSQ_Poly"      = "Bag-STLSQ (Poly)",
    "SINDy_PI_Rational"   = "SINDy-RatPoly2",
    "SINDy_RatPoly2"      = "SINDy-RatPoly2",
    "SR_Operon"           = "SR Pool (MDL)",
    "SR_Pool_MDL"         = "SR Pool (MDL)",
    "SR_PySR"             = "SR Pool (Pareto)",
    "SR_Pool_Pareto"      = "SR Pool (Pareto)",
    "UQ_SINDy_Bayes_Poly" = "BAS-SINDy (Poly)",
    "UQ_SINDy_Bayes_Oracle" = "BAS-SINDy (Oracle)",
    "WSINDy_Poly"         = "WSINDy (Poly)",
    "WSINDy_Oracle"       = "WSINDy (Oracle)"
  )
  out <- ifelse(m_str %in% names(map), map[m_str], m_str)
  unname(out)
}

#' Assign Method Family for Okabe-Ito colorblind-safe palettes
assign_method_family <- function(m_str) {
  m_str <- as.character(m_str)
  ifelse(grepl("^BGNLM", m_str), "BGNLM (Basis-Free)",
    ifelse(grepl("^BGLM", m_str), "Predefined BGLM",
      ifelse(grepl("SINDy|STLSQ|WSINDy|BAS", m_str), "SINDy / UQ-SINDy",
        "Symbolic Regression")))
}

okabe_ito_family_pal <- c(
  "BGNLM (Basis-Free)"  = "#0072B2",
  "Predefined BGLM"     = "#009E73",
  "SINDy / UQ-SINDy"    = "#D55E00",
  "Symbolic Regression" = "#9467BD"
)

#' Clean term label for publication LaTeX tables
clean_term_for_latex <- function(term_str) {
  s <- gsub("^Method", "", as.character(term_str))
  s <- gsub("SystemGroupGroup_", "Group_", s, fixed = TRUE)
  s <- gsub("DiffMethodtvregdiff", "Diff: TVRegDiff", s, fixed = TRUE)
  if (exists("relabel_method_str", mode = "function")) {
    s <- relabel_method_str(s)
  }
  s <- gsub("_", "\\\\_", s)
  s
}

#' Clean metric name for LaTeX tables
clean_metric_latex <- function(m_str) {
  map <- c(
    "F1" = "$F_1$",
    "FDR" = "FDR",
    "TPR" = "TPR",
    "R2_test_ood" = "$R^2_{\\text{OOD}}$",
    "R2_test_in" = "$R^2_{\\text{in}}$",
    "R2_traj" = "$R^2_{\\text{traj}}$",
    "NRMSE_traj" = "$\\text{NRMSE}_{\\text{traj}}$",
    "Exact_Recovery_Logit" = "$I_{\\text{exact}}$ (Logit)",
    "F1_Regime" = "$F_1$ Interaction",
    "Group_A_Polynomial" = "Group A ($F_1$)",
    "Group_B_NonPolynomial" = "Group B ($F_1$)",
    "Group_C_HighDim" = "Group C ($F_1$)"
  )
  out <- ifelse(m_str %in% names(map), map[m_str], gsub("_", "\\\\_", as.character(m_str)))
  unname(out)
}

#' Export Fixed-Effects Mixed-Model Contrasts to Clean Booktabs LaTeX Table
export_mixed_effects_latex_table <- function(fixed_dt, out_path, caption_str, label_str) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(fixed_dt))
  if (!"df" %in% names(dt)) dt[, df := NA_real_]
  if (!"p.value" %in% names(dt)) dt[, p.value := NA_real_]
  if (!"p.dunnett" %in% names(dt)) dt[, p.dunnett := NA_real_]
  if (!"conf.low" %in% names(dt)) dt[, conf.low := estimate - 1.96 * std.error]
  if (!"conf.high" %in% names(dt)) dt[, conf.high := estimate + 1.96 * std.error]

  has_dunnett <- any(!is.na(dt$p.dunnett))
  is_z_model <- all(is.na(dt$df) | !is.finite(dt$df)) || any(grepl("Logit|glmm|ordbeta", dt$Metric, ignore.case = TRUE))
  stat_hdr <- if (is_z_model) "$z$" else "$t$"
  col_spec <- if (has_dunnett) "llccccccc" else "llcccccc"
  hdr_extra <- if (has_dunnett) " & $p_{\\text{adj}}$" else ""

  lines <- c(
    "\\begin{table*}[t]",
    "\\centering",
    paste0("\\caption{", caption_str, "}"),
    paste0("\\label{", label_str, "}"),
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3.5pt}",
    "\\renewcommand{\\arraystretch}{0.88}",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("Metric & Term / Contrast (vs.\\ \\texttt{BGNLM\\_Full\\_Tuned}) & $\\hat{\\beta}$ & SE & $df$ & ", stat_hdr, " & $95\\%$ CI & $p$-value", hdr_extra, " \\\\"),
    "\\midrule"
  )

  prev_m <- ""
  for (i in seq_len(nrow(dt))) {
    r <- dt[i]
    cur_m <- as.character(r$Metric)
    if (i > 1L && cur_m != prev_m) {
      lines <- c(lines, "\\midrule")
    }
    prev_m <- cur_m

    term_clean <- clean_term_for_latex(r$term)
    metric_clean <- clean_metric_latex(cur_m)
    df_str <- if (is.na(r$df) || !is.finite(r$df)) "--" else sprintf("%.1f", r$df)
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

  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table*}")
  writeLines(lines, out_path)
  invisible(out_path)
}

#' Export ICC Variance Decomposition Table to LaTeX (Booktabs, no vertical rules)
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
    "\\setlength{\\tabcolsep}{3.5pt}",
    paste0("\\begin{tabular}{", col_spec, "}"),
    "\\toprule",
    paste0("Model & Random Effect Group & $\\hat{\\sigma}^2$ & $\\hat{\\sigma}$ & ICC", sing_hdr, " \\\\"),
    "\\midrule"
  )
  prev_mod <- ""
  for (i in seq_len(nrow(dt))) {
    r <- dt[i]
    if (i > 1L && r$Model != prev_mod) lines <- c(lines, "\\midrule")
    prev_mod <- r$Model
    m_clean <- gsub("_", "\\\\_", r$Model)
    g_clean <- gsub("_", "\\\\_", r$Group)
    s_col <- if (has_sing) paste0(" & ", r$Singular) else ""
    lines <- c(lines, sprintf(
      "%s & \\texttt{%s} & %.4f & %.4f & %.3f%s \\\\",
      m_clean, g_clean, r$Variance, r$StdDev, r$ICC, s_col
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, out_path)
  invisible(out_path)
}

#' Generate Publication-Quality Mixed-Effects Forest Plot (Dual Unadjusted & Simultaneous Dunnett 95% CIs)
plot_mixed_effects_forest <- function(fixed_dt, out_pdf, title_str = NULL, dunnett_dt = NULL) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(fixed_dt))
  dt <- dt[grepl("^Method", term)]
  if (nrow(dt) == 0L) return(invisible(NULL))

  dt[, Raw_Method := gsub("^Method", "", term)]
  dt[, Method_Label := human_method_label(Raw_Method)]
  dt[, Family := factor(assign_method_family(Raw_Method), levels = names(okabe_ito_family_pal))]

  if (!"conf.low" %in% names(dt)) dt[, conf.low := estimate - 1.96 * std.error]
  if (!"conf.high" %in% names(dt)) dt[, conf.high := estimate + 1.96 * std.error]

  dt[, dun.low := conf.low]
  dt[, dun.high := conf.high]
  if (!is.null(dunnett_dt) && nrow(dunnett_dt) > 0L && all(c("dunnett.low", "dunnett.high") %in% names(dunnett_dt))) {
    dd <- copy(as.data.table(dunnett_dt))
    dd[, Raw_Method := sub(" - .*", "", contrast)]
    for (i in seq_len(nrow(dt))) {
      hit <- dd[Metric == dt$Metric[i] & Raw_Method == dt$Raw_Method[i]]
      if (nrow(hit) >= 1L) {
        dt[i, dun.low := hit$dunnett.low[1L]]
        dt[i, dun.high := hit$dunnett.high[1L]]
        if ("unadj.low" %in% names(hit) && is.finite(hit$unadj.low[1L])) {
          dt[i, conf.low := hit$unadj.low[1L]]
          dt[i, conf.high := hit$unadj.high[1L]]
        }
      }
    }
  }

  metric_map <- c(
    "F1" = "Structural F1 Contrast",
    "FDR" = "False Discovery Rate Contrast",
    "R2_test_ood" = "Derivative R2 (OOD) Contrast",
    "R2_traj" = "Trajectory Rollout R2 Contrast",
    "NRMSE_traj" = "Rollout NRMSE Contrast"
  )
  ordered_metrics <- unique(c(intersect(names(metric_map), unique(dt$Metric)),
                              setdiff(unique(dt$Metric), names(metric_map))))
  panel_labels <- setNames(
    paste0("(", letters[seq_along(ordered_metrics)], ") ",
           ifelse(ordered_metrics %in% names(metric_map), metric_map[ordered_metrics], ordered_metrics)),
    ordered_metrics
  )
  dt[, Metric_Panel := factor(panel_labels[Metric], levels = unname(panel_labels))]

  # Order methods consistently by F1 (or first metric) contrast
  ord_df <- dt[Metric == dt$Metric[1], .(ord_est = mean(estimate)), by = Method_Label]
  dt <- merge(dt, ord_df, by = "Method_Label", all.x = TRUE)
  dt[, Method_Label := reorder(Method_Label, ord_est)]

  n_panels <- uniqueN(dt$Metric_Panel)
  nrow_wrap <- if (n_panels == 4L) 2L else 1L
  fig_w <- 7.4
  fig_h <- if (n_panels == 4L) 4.6 else 3.0

  p <- ggplot(dt, aes(x = estimate, y = Method_Label, color = Family)) +
    geom_vline(xintercept = 0.0, linetype = "dashed", color = "grey30", linewidth = 0.55) +
    geom_errorbar(aes(xmin = dun.low, xmax = dun.high), orientation = "y", width = 0.30, linewidth = 0.50, alpha = 0.85) +
    geom_segment(aes(x = conf.low, xend = conf.high, y = Method_Label, yend = Method_Label), linewidth = 1.25) +
    geom_point(size = 2.0) +
    scale_color_manual(values = okabe_ito_family_pal, drop = TRUE) +
    facet_wrap(~ Metric_Panel, nrow = nrow_wrap, scales = "free_x") +
    labs(
      x = "Mixed-Effects Contrast vs. BGNLM (Tuned) [Inner: Unadj. 95% CI, Outer: Simultaneous Dunnett 95% CI]",
      y = NULL,
      color = "Method Family"
    ) +
    theme_bw(base_size = 8.5) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 7.8, face = "bold"),
      legend.text = element_text(size = 7.0),
      strip.background = element_rect(fill = "#eef2f6", color = "grey75"),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7.5),
      axis.title.x = element_text(size = 7.6),
      plot.margin = margin(5, 10, 5, 5)
    ) +
    guides(color = guide_legend(nrow = 1))

  ggsave(out_pdf, p, width = fig_w, height = fig_h, bg = "white")
  ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = fig_w, height = fig_h, dpi = 220, bg = "white")
  paper_fig <- file.path("paper_icml/figures", basename(out_pdf))
  if (dir.exists("paper_icml/figures")) ggsave(paper_fig, p, width = fig_w, height = fig_h, bg = "white")
  invisible(p)
}

#' Generate 4-Panel Full-Width Benchmark Overview (F1, FDR, R2_OOD, R2_traj across Groups A/B/C)
plot_benchmark_summary <- function(df_synth, out_pdf = "results/figures/fig_synth_overview.pdf") {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  dt <- copy(as.data.table(df_synth))
  if (nrow(dt) == 0L) return(invisible(NULL))

  dt[, Method_Label := human_method_label(Method)]
  group_map <- c(
    "Group_A_Polynomial"    = "Group A: Poly",
    "Group_B_NonPolynomial" = "Group B: Non-Poly",
    "Group_C_HighDim"       = "Group C: 5D Lorenz-96"
  )
  dt[, Group_Short := factor(group_map[SystemGroup], levels = unname(group_map))]

  # Aggregate across all reps per Method x Group
  long_dt <- rbindlist(list(
    dt[, .(Val = mean(F1, na.rm = TRUE), SE = sd(F1, na.rm = TRUE) / sqrt(sum(!is.na(F1))), Panel = "(a) Structural F1 (higher better)"), by = .(Method_Label, Group_Short)],
    dt[, .(Val = mean(FDR, na.rm = TRUE), SE = sd(FDR, na.rm = TRUE) / sqrt(sum(!is.na(FDR))), Panel = "(b) False Discovery Rate (lower better)"), by = .(Method_Label, Group_Short)],
    dt[, .(Val = mean(R2_test_ood, na.rm = TRUE), SE = sd(R2_test_ood, na.rm = TRUE) / sqrt(sum(!is.na(R2_test_ood))), Panel = "(c) Derivative R2 OOD (higher better)"), by = .(Method_Label, Group_Short)],
    dt[, .(Val = mean(R2_traj, na.rm = TRUE), SE = sd(R2_traj, na.rm = TRUE) / sqrt(sum(!is.na(R2_traj))), Panel = "(d) Trajectory Rollout R2 (higher better)"), by = .(Method_Label, Group_Short)]
  ))
  long_dt[, Panel := factor(Panel, levels = c(
    "(a) Structural F1 (higher better)",
    "(b) False Discovery Rate (lower better)",
    "(c) Derivative R2 OOD (higher better)",
    "(d) Trajectory Rollout R2 (higher better)"
  ))]

  ord_m <- long_dt[Panel == "(a) Structural F1 (higher better)", .(m = mean(Val, na.rm = TRUE)), by = Method_Label]
  long_dt <- merge(long_dt, ord_m, by = "Method_Label", all.x = TRUE)
  long_dt[, Method_Label := reorder(Method_Label, m)]

  grp_pal <- c(
    "Group A: Poly" = "#0072B2",
    "Group B: Non-Poly" = "#009E73",
    "Group C: 5D Lorenz-96" = "#D55E00"
  )

  p <- ggplot(long_dt, aes(x = Method_Label, y = Val, fill = Group_Short)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    geom_errorbar(
      aes(ymin = pmax(0, Val - 1.96 * SE), ymax = pmin(1, Val + 1.96 * SE)),
      position = position_dodge(width = 0.78), width = 0.25, linewidth = 0.35
    ) +
    coord_flip(ylim = c(0, 1.02)) +
    scale_fill_manual(values = grp_pal) +
    facet_wrap(~ Panel, nrow = 2, ncol = 2) +
    labs(x = NULL, y = "Mean Score (+/- 1.96 SE)", fill = "Dynamical Regime") +
    theme_bw(base_size = 8.5) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7.5),
      strip.background = element_rect(fill = "#eef2f6", color = "grey75"),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 7.2)
    )

  ggsave(out_pdf, p, width = 6.75, height = 4.6, bg = "white")
  ggsave(sub("\\.pdf$", ".png", out_pdf), p, width = 6.75, height = 4.6, dpi = 220, bg = "white")
  if (dir.exists("paper_icml/figures")) ggsave(file.path("paper_icml/figures", basename(out_pdf)), p, width = 6.75, height = 4.6, bg = "white")
  invisible(p)
}
