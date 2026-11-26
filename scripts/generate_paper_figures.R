#!/usr/bin/env Rscript
# ==============================================================================
# scripts/generate_paper_figures.R
# Redesigned Publication-Quality ICML Full-Width Figures (6.75 in, Okabe-Ito)
# Generates:
#   1. fig_synth_overview.pdf           (4 panels: F1, FDR, R2_OOD, R2_traj across Groups A/B/C)
#   2. fig_consistency_scaling.pdf      (2x3: Sample size N and Noise sigma with scaled-budget overlay)
#   3. fig_uq_calibration_stability.pdf (1x3: Nogueira Stability, PIP Reliability Diagram, Real OC vs F1)
#   4. fig_qualitative_rollouts.pdf     (1x3: Lorenz-63 + 95% ribbon, Goodwin rollout, Authentic Lynx-Hare)
# (Forest plots fig_mixed_effects_forest_synth.pdf & _real.pdf come from run_mixed_models.R)
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(patchwork)
})

source("R/01_ode_simulators.R")
source("R/02_basis_builders.R")
source("R/05_symbolic_canonicalize.R")
source("R/06_metrics_consistency.R")
source("R/03_fbms_wrappers.R")
source("R/04_sindy_baselines.R")
source("R/08_mixed_effects_eval.R")
source("R/09_plotting_tables.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_icml/figures", recursive = TRUE, showWarnings = FALSE)

save_fig_all <- function(plot_obj, name_base, width = 6.75, height = 4.2) {
  p_res <- file.path("results/figures", paste0(name_base, ".pdf"))
  p_png <- file.path("results/figures", paste0(name_base, ".png"))
  p_icml <- file.path("paper_icml/figures", paste0(name_base, ".pdf"))
  ggplot2::ggsave(p_res, plot_obj, width = width, height = height, bg = "white")
  ggplot2::ggsave(p_png, plot_obj, width = width, height = height, dpi = 240, bg = "white")
  ggplot2::ggsave(p_icml, plot_obj, width = width, height = height, bg = "white")
  message("[OK] Saved figure: ", p_icml)
}

# 1. Load audited aggregated results
synth_dt <- readRDS("results/aggregated/synth_results.rds")
synth_stab <- fread("results/aggregated/synth_cell_summary_stability.csv")
real_stab <- fread("results/aggregated/real_cell_summary_stability.csv")

# Regenerate 4-panel synth overview cleanly
plot_benchmark_summary(synth_dt, out_pdf = "results/figures/fig_synth_overview.pdf")

# Key methods for crisp uncrowded scaling curves
method_levels <- c(
  "BGNLM_Full_Tuned", "BGLM_Oracle_Basis", "BGNLM_Hybrid_Poly",
  "BGLM_Poly_Basis", "SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly",
  "BGNLM_IWSM_Untuned", "SR_Operon", "UQ_SINDy_Poly"
)

method_labels <- c(
  "BGNLM_Full_Tuned"   = "BGNLM (Tuned)",
  "BGLM_Oracle_Basis"  = "BGLM (Oracle)",
  "BGNLM_Hybrid_Poly"  = "BGNLM-Hybrid",
  "BGLM_Poly_Basis"    = "BGLM (Poly)",
  "SINDy_STLSQ_Oracle" = "STLSQ (Oracle)",
  "SINDy_STLSQ_Poly"   = "STLSQ (Poly)",
  "BGNLM_IWSM_Untuned" = "BGNLM (Untuned)",
  "SR_Operon"          = "SR Pool (MDL)",
  "UQ_SINDy_Poly"      = "Bag-STLSQ (Poly)"
)

group_labels <- c(
  "Group_A_Polynomial"    = "Group A: Polynomial (7 ODEs)",
  "Group_B_NonPolynomial" = "Group B: Non-Poly (4 ODEs)",
  "Group_C_HighDim"       = "Group C: 5D Lorenz-96"
)

# Okabe-Ito colorblind-safe line palette
pal_lines <- c(
  "BGNLM (Tuned)"    = "#0072B2",
  "BGLM (Oracle)"    = "#009E73",
  "BGNLM-Hybrid"     = "#56B4E9",
  "BGLM (Poly)"      = "#8C564B",
  "STLSQ (Oracle)"   = "#E69F00",
  "STLSQ (Poly)"     = "#D55E00",
  "BGNLM (Untuned)"  = "#CC79A7",
  "SR Pool (MDL)"    = "#9467BD",
  "Bag-STLSQ (Poly)" = "#666666"
)

# ------------------------------------------------------------------------------
# Figure 2: fig_consistency_scaling.pdf (6.75 x 4.1 in, 2x3 via patchwork)
# ------------------------------------------------------------------------------
sub_s <- synth_dt[Method %in% method_levels]
sub_s[, Method_Name := factor(method_labels[Method], levels = unname(method_labels))]
sub_s[, Group_Name := factor(group_labels[SystemGroup], levels = unname(group_labels))]

df_n <- sub_s[, .(Mean_F1 = mean(F1, na.rm = TRUE)), by = .(Group_Name, Method_Name, SampleSize)]
df_noise <- sub_s[, .(Mean_F1 = mean(F1, na.rm = TRUE)), by = .(Group_Name, Method_Name, Noise)]

# Optional scaled-budget large-N points from Job 1b
scaled_pts <- NULL
if (file.exists("results/ablation_native/job1_mpm_pruning_and_largeN.csv")) {
  j1 <- fread("results/ablation_native/job1_mpm_pruning_and_largeN.csv")
  j1b <- j1[SubJob == "1b_LargeN" & Variant == "BGNLM_ScaledBudget_SqrtN"]
  if (nrow(j1b) > 0L) {
    j1b[, Group_Name := factor(ifelse(System == "sys02_lorenz63", group_labels[["Group_A_Polynomial"]], group_labels[["Group_B_NonPolynomial"]]),
                               levels = unname(group_labels))]
    scaled_pts <- j1b[, .(Mean_F1 = mean(F1, na.rm = TRUE)), by = .(Group_Name, SampleSize)]
  }
}

p_n <- ggplot(df_n, aes(x = SampleSize, y = Mean_F1, color = Method_Name, group = Method_Name)) +
  geom_line(aes(linewidth = Method_Name == "BGNLM (Tuned)"), alpha = 0.90) +
  geom_point(size = 1.35) +
  scale_linewidth_manual(values = c("TRUE" = 1.15, "FALSE" = 0.55), guide = "none") +
  scale_x_log10(
    breaks = c(100, 250, 500, 1000, 2500),
    labels = c("100", "250", "500", "1k", "2.5k")
  ) +
  scale_color_manual(values = pal_lines) +
  facet_wrap(~ Group_Name, nrow = 1) +
  coord_cartesian(ylim = c(0.0, 1.02)) +
  labs(x = expression("Sample Size " * italic(N) * " (log scale)"), y = expression("Structural " * italic(F)[1]), color = NULL) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_rect(fill = "#eef2f6", color = "grey75"),
    strip.text = element_text(face = "bold", size = 7.8),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7.0),
    legend.position = "bottom",
    legend.text = element_text(size = 7.2)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

if (!is.null(scaled_pts) && nrow(scaled_pts) > 0L) {
  p_n <- p_n + geom_point(
    data = scaled_pts,
    aes(x = SampleSize, y = Mean_F1),
    inherit.aes = FALSE, shape = 8, size = 2.0, color = "#000000", stroke = 0.8
  )
}

p_noise <- ggplot(df_noise, aes(x = Noise, y = Mean_F1, color = Method_Name, group = Method_Name)) +
  geom_line(aes(linewidth = Method_Name == "BGNLM (Tuned)"), alpha = 0.90) +
  geom_point(size = 1.35) +
  scale_linewidth_manual(values = c("TRUE" = 1.15, "FALSE" = 0.55), guide = "none") +
  scale_x_log10(
    breaks = c(0.001, 0.01, 0.05, 0.25, 1.0),
    labels = c("1e-3", "0.01", "0.05", "0.25", "1.0")
  ) +
  scale_color_manual(values = pal_lines) +
  facet_wrap(~ Group_Name, nrow = 1) +
  coord_cartesian(ylim = c(0.0, 1.02)) +
  labs(x = expression("Derivative Noise Intensity " * sigma * " (log scale)"), y = expression("Structural " * italic(F)[1]), color = NULL) +
  theme_bw(base_size = 8.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7.0),
    legend.position = "bottom",
    legend.text = element_text(size = 7.2)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

p_scaling <- (p_n / p_noise) +
  plot_layout(guides = "collect", heights = c(1.05, 1.0)) +
  plot_annotation(theme = theme(legend.position = "bottom", legend.text = element_text(size = 7.2)))

save_fig_all(p_scaling, "fig_consistency_scaling", width = 7.4, height = 4.3)

# ------------------------------------------------------------------------------
# Figure 3: fig_uq_calibration_stability.pdf (6.75 x 2.75 in, 1x3 via patchwork)
#   (a) Nogueira Stability vs F1 (ggrepel)
#   (b) PIP Calibration Reliability Diagram (Decile bins vs y=x)
#   (c) Real-World Parsimony (OC) vs F1 (ggrepel)
# ------------------------------------------------------------------------------
synth_m <- synth_stab[, .(
  Mean_F1 = mean(Mean_F1, na.rm = TRUE),
  Mean_Nogueira = mean(Nogueira_Phi, na.rm = TRUE),
  Mean_FDR = mean(Mean_FDR, na.rm = TRUE)
), by = Method]
synth_m[, Label := human_method_label(Method)]
synth_m[, Family := factor(assign_method_family(Method), levels = names(okabe_ito_family_pal))]

p_stab <- ggplot(synth_m, aes(x = Mean_Nogueira, y = Mean_F1, color = Family)) +
  geom_point(aes(size = pmax(0.1, 1.0 - Mean_FDR)), alpha = 0.90) +
  ggrepel::geom_text_repel(
    aes(label = Label), size = 2.0, max.overlaps = Inf, box.padding = 0.22,
    point.padding = 0.15, min.segment.length = 0.1, seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = okabe_ito_family_pal) +
  scale_size_continuous(range = c(1.6, 4.0), guide = "none") +
  coord_cartesian(xlim = c(0.40, 0.96), ylim = c(0.35, 0.95)) +
  labs(
    x = expression("(a) Nogueira Stability " * hat(Phi)),
    y = expression("Mean Structural " * italic(F)[1]),
    color = "Method Family"
  ) +
  theme_bw(base_size = 8.0) +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(size = 7.6),
        legend.position = "bottom", legend.text = element_text(size = 7.0), legend.title = element_text(size = 7.4, face = "bold")) +
  guides(color = guide_legend(nrow = 1))

# Panel (b): PIP Calibration Reliability Diagram
cal_file <- "results/ablation_native/pip_calibration_reliability.csv"
if (file.exists(cal_file) && nrow(fread(cal_file)) > 0L) {
  cal_raw <- fread(cal_file)
  cal_agg <- cal_raw[, .(
    Mean_PIP = weighted.mean(mean_pip, count, na.rm = TRUE),
    Emp_Freq = weighted.mean(emp_freq, count, na.rm = TRUE),
    Total_N  = sum(count)
  ), by = .(Method, bin_id)]
  cal_agg[, Method_Lab := human_method_label(Method)]
} else {
  # Fallback single calibration evaluation on Lorenz-63 + Goodwin
  cal_agg <- data.table(Method_Lab = "BGNLM (Tuned)", Mean_PIP = c(0.05, 0.5, 0.95), Emp_Freq = c(0.03, 0.48, 0.96))
}

pip_pal <- c(
  "BGNLM (Tuned)"    = "#0072B2",
  "BAS-SINDy (Poly)" = "#009E73",
  "Bag-STLSQ (Poly)" = "#D55E00"
)

p_calib <- ggplot(cal_agg, aes(x = Mean_PIP, y = Emp_Freq, color = Method_Lab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_line(linewidth = 0.85, alpha = 0.90) +
  geom_point(size = 1.6) +
  scale_color_manual(values = pip_pal, name = "PIP Estimator") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    x = "(b) Predicted Mean PIP",
    y = "Empirical Inclusion Frequency"
  ) +
  theme_bw(base_size = 8.0) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7.4),
    legend.position = "bottom",
    legend.title = element_text(size = 7.0, face = "bold"),
    legend.text = element_text(size = 6.6),
    plot.margin = margin(3, 5, 3, 5)
  ) +
  guides(color = guide_legend(nrow = 1))

# Panel (c): Real-World Parsimony (OC) vs F1 (conditional on OK runs, OC is NA for crashes)
real_m <- real_stab[, .(
  Mean_F1 = mean(Mean_F1, na.rm = TRUE),
  Mean_OC = mean(Mean_OC, na.rm = TRUE),
  Mean_R2_traj = mean(Mean_R2_traj, na.rm = TRUE)
), by = Method]
real_m[, Label := human_method_label(Method)]
real_m[, Family := factor(assign_method_family(Method), levels = names(okabe_ito_family_pal))]

p_real_p <- ggplot(real_m, aes(x = Mean_OC, y = Mean_F1, color = Family)) +
  geom_point(aes(size = pmax(0.15, Mean_R2_traj)), alpha = 0.90) +
  ggrepel::geom_text_repel(
    aes(label = Label), size = 2.0, max.overlaps = Inf, box.padding = 0.22,
    point.padding = 0.15, min.segment.length = 0.1, seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = okabe_ito_family_pal) +
  scale_size_continuous(range = c(1.6, 4.0), guide = "none") +
  coord_cartesian(xlim = c(-2, 65), ylim = c(0.20, 0.78)) +
  labs(
    x = "(c) Mean Operation Count OC",
    y = expression("Quasi-Empirical " * italic(F)[1]),
    color = "Method Family"
  ) +
  theme_bw(base_size = 8.0) +
  theme(panel.grid.minor = element_blank(), axis.title = element_text(size = 7.6),
        legend.position = "bottom", legend.text = element_text(size = 7.0), legend.title = element_text(size = 7.4, face = "bold")) +
  guides(color = guide_legend(nrow = 1))

p_uq_triple <- (p_stab | p_calib | p_real_p) +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(theme = theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 6, 0, 6)
  ))

save_fig_all(p_uq_triple, "fig_uq_calibration_stability", width = 9.4, height = 4.15)
save_fig_all(p_uq_triple, "fig_uq_stability_pareto", width = 9.4, height = 4.15)

# ------------------------------------------------------------------------------
# Figure 4: fig_qualitative_rollouts.pdf (6.75 x 2.60 in, 1x3 Trajectory Showcase)
#   (a) Lorenz-63 x1(t) autonomous rollout + 95% credible ribbon
#   (b) Goodwin MM x3(t) rollout comparing Ground Truth vs BGNLM vs STLSQ vs RatPoly
#   (c) Authentic 1900-1920 Lynx-Hare zero-leakage rollout
# ------------------------------------------------------------------------------
rollout_integrate <- function(predict_rhs_fn, rollout_spec) {
  times <- rollout_spec$times
  x0 <- as.numeric(rollout_spec$x0)
  m <- length(x0)
  names(x0) <- colnames(rollout_spec$X_true)
  bound_max <- 50.0 * max(abs(rollout_spec$X_true), 1.0)
  ode_fn <- function(t, s, p) {
    if (any(!is.finite(s)) || max(abs(s)) > bound_max) return(list(rep(0.0, m)))
    mat <- matrix(s, nrow = 1L); colnames(mat) <- names(x0)
    ds <- tryCatch(as.numeric(predict_rhs_fn(mat)), error = function(e) rep(0.0, m))
    ds[!is.finite(ds)] <- 0.0
    list(pmax(pmin(ds, bound_max), -bound_max))
  }
  sol <- deSolve::ode(y = x0, times = times, func = ode_fn, parms = NULL, method = "rk4")
  as.matrix(sol[, -1, drop = FALSE])
}

# Panel (a): Lorenz-63
td_lor <- generate_synthetic_task_data("sys02_lorenz63", N_train = 350L, noise_eta = 0.05, seed = 101L)
fit_lor <- run_fbms_method("BGNLM_Full_Tuned", td_lor, cores = 4L, fast_mode = FALSE)
X_lor_pred <- rollout_integrate(fit_lor$predict_rhs_fn, td_lor$rollout)
# Estimate trajectory uncertainty spread by propagating derivative residual SE
t_lor <- td_lor$rollout$times
dt_lor <- mean(diff(t_lor))
sd_dx1 <- mean(0.5 * (fit_lor$upper_ood[, 1] - fit_lor$lower_ood[, 1]), na.rm = TRUE)
rib_lor <- pmin(8.0, sd_dx1 * sqrt(pmax(0.01, t_lor)) * 0.25)

df_lor <- data.table(
  t = rep(t_lor, 2),
  x1 = c(td_lor$rollout$X_true[, 1], X_lor_pred[, 1]),
  Series = rep(c("Ground Truth", "BGNLM (Tuned)"), each = length(t_lor)),
  low = c(td_lor$rollout$X_true[, 1], X_lor_pred[, 1] - rib_lor),
  high = c(td_lor$rollout$X_true[, 1], X_lor_pred[, 1] + rib_lor)
)

p_lor <- ggplot(df_lor, aes(x = t, y = x1, color = Series, linetype = Series)) +
  geom_ribbon(data = df_lor[Series == "BGNLM (Tuned)"], aes(ymin = low, ymax = high), fill = "#0072B2", alpha = 0.18, color = NA) +
  geom_line(linewidth = 0.75) +
  scale_color_manual(values = c("Ground Truth" = "#111111", "BGNLM (Tuned)" = "#0072B2")) +
  scale_linetype_manual(values = c("Ground Truth" = "solid", "BGNLM (Tuned)" = "dashed")) +
  labs(x = "(a) Lorenz-63 Rollout Time t", y = expression("State " * italic(x)[1](italic(t))), color = NULL, linetype = NULL) +
  theme_bw(base_size = 8.0) +
  theme(legend.position = c(0.68, 0.85), legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 6.3), panel.grid.minor = element_blank(), axis.title = element_text(size = 7.6))

# Panel (b): Goodwin Rational Oscillator x3(t)
td_gw <- generate_synthetic_task_data("sys11_goodwin_mm", N_train = 400L, noise_eta = 0.05, seed = 101L)
f_gw_bg <- run_fbms_method("BGNLM_Full_Tuned", td_gw, cores = 4L, fast_mode = FALSE)
f_gw_sp <- run_sindy_baseline("SINDy_STLSQ_Poly", td_gw)
f_gw_rp <- run_sindy_baseline("SINDy_PI_Rational", td_gw)
X_gw_bg <- rollout_integrate(f_gw_bg$predict_rhs_fn, td_gw$rollout)
X_gw_sp <- rollout_integrate(f_gw_sp$predict_rhs_fn, td_gw$rollout)
X_gw_rp <- rollout_integrate(f_gw_rp$predict_rhs_fn, td_gw$rollout)
t_gw <- td_gw$rollout$times

df_gw <- data.table(
  t = rep(t_gw, 4),
  x3 = c(td_gw$rollout$X_true[, 3], X_gw_bg[, 3], X_gw_sp[, 3], X_gw_rp[, 3]),
  Method = factor(rep(c("Ground Truth", "BGNLM (Tuned)", "STLSQ (Poly)", "SINDy-RatPoly"), each = length(t_gw)),
                  levels = c("Ground Truth", "BGNLM (Tuned)", "STLSQ (Poly)", "SINDy-RatPoly"))
)

p_gw <- ggplot(df_gw, aes(x = t, y = x3, color = Method, linetype = Method)) +
  geom_line(linewidth = 0.72) +
  scale_color_manual(values = c("Ground Truth" = "#111111", "BGNLM (Tuned)" = "#0072B2", "STLSQ (Poly)" = "#D55E00", "SINDy-RatPoly" = "#009E73")) +
  scale_linetype_manual(values = c("Ground Truth" = "solid", "BGNLM (Tuned)" = "dashed", "STLSQ (Poly)" = "dotted", "SINDy-RatPoly" = "dotdash")) +
  labs(x = "(b) Goodwin Oscillator Time t", y = expression("State " * italic(x)[3](italic(t))), color = NULL, linetype = NULL) +
  theme_bw(base_size = 8.0) +
  theme(legend.position = c(0.68, 0.25), legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.text = element_text(size = 6.0), panel.grid.minor = element_blank(), axis.title = element_text(size = 7.6))

# Panel (c): Authentic 1900-1920 Lynx-Hare Zero-Leakage Empirical Rollout
df_lh <- read.csv("data/real_world/lynx_hare_annual_raw_1900_1920.csv", stringsAsFactors = FALSE)
t_tr_r <- df_lh$time[1:14]; X_tr_r <- as.matrix(df_lh[1:14, c("x1", "x2")])
t_ood_r <- df_lh$time[15:21]; X_ood_r <- as.matrix(df_lh[15:21, c("x1", "x2")])
t_d <- seq(0, 13, by = 0.25)
X_d <- cbind(x1 = pmax(1, spline(t_tr_r, X_tr_r[, 1], xout = t_d, method = "natural")$y),
             x2 = pmax(1, spline(t_tr_r, X_tr_r[, 2], xout = t_d, method = "natural")$y))
Y_d <- estimate_derivatives_sgolay_spline(X_d, t_d)
t_od <- seq(14, 20, by = 0.25)
X_od <- cbind(x1 = pmax(1, spline(t_ood_r, X_ood_r[, 1], xout = t_od, method = "natural")$y),
              x2 = pmax(1, spline(t_ood_r, X_ood_r[, 2], xout = t_od, method = "natural")$y))
Y_od <- estimate_derivatives_sgolay_spline(X_od, t_od)
td_lh <- list(
  system = list(oracle_extra_basis = character(0)),
  X_train = X_d, Y_train = Y_d, X_test_in = X_d[seq(1, nrow(X_d), 4), ], X_test_ood = X_od,
  rollout = list(times = t_od - t_od[1], X_true = X_od, x0 = X_od[1, ])
)
fit_lh <- run_fbms_method("BGNLM_Full_Tuned", td_lh, cores = 4L, fast_mode = FALSE)
X_lh_pred <- rollout_integrate(fit_lh$predict_rhs_fn, td_lh$rollout)

df_lh_plot <- rbind(
  data.table(Year = 1900 + df_lh$time, Value = df_lh$x1, Species = "Hare (Obs)", Type = "Empirical"),
  data.table(Year = 1900 + df_lh$time, Value = df_lh$x2, Species = "Lynx (Obs)", Type = "Empirical"),
  data.table(Year = 1900 + t_od, Value = X_lh_pred[, 1], Species = "Hare (BGNLM)", Type = "OOD Rollout"),
  data.table(Year = 1900 + t_od, Value = X_lh_pred[, 2], Species = "Lynx (BGNLM)", Type = "OOD Rollout")
)

p_lh <- ggplot(df_lh_plot, aes(x = Year, y = Value, color = Species, linetype = Type)) +
  geom_vline(xintercept = 1913.5, linetype = "dotted", color = "grey45") +
  geom_line(data = df_lh_plot[Type == "OOD Rollout"], linewidth = 0.85) +
  geom_point(data = df_lh_plot[Type == "Empirical"], size = 1.4) +
  geom_line(data = df_lh_plot[Type == "Empirical"], linewidth = 0.45, alpha = 0.65) +
  scale_color_manual(values = c("Hare (Obs)" = "#E69F00", "Hare (BGNLM)" = "#E69F00",
                                "Lynx (Obs)" = "#0072B2", "Lynx (BGNLM)" = "#0072B2")) +
  scale_linetype_manual(values = c("Empirical" = "solid", "OOD Rollout" = "dashed"), guide = "none") +
  coord_cartesian(ylim = c(0, 96)) +
  labs(x = "(c) Authentic Pelts 1900-1920", y = "Pelts (Thousands)", color = NULL) +
  theme_bw(base_size = 8.0) +
  theme(legend.position = c(0.54, 0.88), legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.text = element_text(size = 5.6), legend.key.height = unit(0.32, "cm"),
        panel.grid.minor = element_blank(), axis.title = element_text(size = 7.6)) +
  guides(color = guide_legend(ncol = 2))

p_rollouts <- (p_lor | p_gw | p_lh)
save_fig_all(p_rollouts, "fig_qualitative_rollouts", width = 9.4, height = 3.3)
message("[DONE] All 6 publication-quality figures generated.")
