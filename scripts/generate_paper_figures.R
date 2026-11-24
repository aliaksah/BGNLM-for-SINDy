#!/usr/bin/env Rscript
# ==============================================================================
# scripts/generate_paper_figures.R
# Generates Publication-Quality Multi-Panel PDF/PNG Figures for paper_icml/
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("paper_icml/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Load aggregated data
synth_dt <- readRDS("results/aggregated/synth_results.rds")
synth_stab <- fread("results/aggregated/synth_cell_summary_stability.csv")
real_stab <- fread("results/aggregated/real_cell_summary_stability.csv")

# Clean human-readable labels and cohesive color palette
method_levels <- c(
  "BGNLM_Full_Tuned", "BGLM_Oracle_Basis", "BGNLM_Hybrid_Poly",
  "BGLM_Poly_Basis", "SINDy_STLSQ_Oracle", "BGLM_Overcomplete",
  "SINDy_STLSQ_Poly", "BGNLM_IWSM_Untuned", "SR_Operon", "SR_PySR"
)

method_labels <- c(
  "BGNLM_Full_Tuned"   = "BGNLM (Tuned)",
  "BGLM_Oracle_Basis"  = "BGLM (Oracle)",
  "BGNLM_Hybrid_Poly"  = "BGNLM-Hybrid",
  "BGLM_Poly_Basis"    = "BGLM (Poly)",
  "SINDy_STLSQ_Oracle" = "STLSQ (Oracle)",
  "BGLM_Overcomplete"  = "BGLM (Overcomplete)",
  "SINDy_STLSQ_Poly"   = "STLSQ (Poly)",
  "BGNLM_IWSM_Untuned" = "BGNLM (Untuned IWSM)",
  "SR_Operon"          = "SR Pool (MDL)",
  "SR_PySR"            = "SR Pool (Pareto)"
)

group_labels <- c(
  "Group_A_Polynomial"    = "Group A: Polynomial ODEs (7 Systems)",
  "Group_B_NonPolynomial" = "Group B: Non-Polynomial ODEs (4 Systems)",
  "Group_C_HighDim"       = "Group C: High-Dim 5D Lorenz-96"
)

pal <- c(
  "BGNLM (Tuned)"        = "#0072B2",
  "BGLM (Oracle)"        = "#009E73",
  "BGNLM-Hybrid"         = "#56B4E9",
  "BGLM (Poly)"          = "#8C564B",
  "STLSQ (Oracle)"       = "#E69F00",
  "BGLM (Overcomplete)"  = "#7F7F7F",
  "STLSQ (Poly)"         = "#D55E00",
  "BGNLM (Untuned IWSM)" = "#CC79A7",
  "SR Pool (MDL)"        = "#9467BD",
  "SR Pool (Pareto)"     = "#BCBD22"
)

# ------------------------------------------------------------------------------
# Figure A: fig_consistency_scaling.pdf (3-panel Sample Size & Noise Scaling)
# ------------------------------------------------------------------------------
sub_s <- synth_dt[Method %in% method_levels]
sub_s[, Method_Name := factor(method_labels[Method], levels = unname(method_labels))]
sub_s[, Group_Name := factor(group_labels[SystemGroup], levels = unname(group_labels))]

df_n <- sub_s[, .(
  Mean_F1 = mean(F1, na.rm = TRUE),
  SE_F1 = sd(F1, na.rm = TRUE) / sqrt(.N)
), by = .(Group_Name, Method_Name, SampleSize)]
df_n[, Axis_Type := "Sample Size N (log scale)"]
df_n[, X_val := as.numeric(SampleSize)]

df_noise <- sub_s[, .(
  Mean_F1 = mean(F1, na.rm = TRUE),
  SE_F1 = sd(F1, na.rm = TRUE) / sqrt(.N)
), by = .(Group_Name, Method_Name, Noise)]
df_noise[, Axis_Type := "Noise Level sigma (log scale)"]
df_noise[, X_val := as.numeric(Noise)]

p_n <- ggplot(df_n, aes(x = SampleSize, y = Mean_F1, color = Method_Name, group = Method_Name)) +
  geom_line(aes(linewidth = Method_Name == "BGNLM (Tuned)"), alpha = 0.92) +
  geom_point(size = 1.8) +
  scale_linewidth_manual(values = c("TRUE" = 1.35, "FALSE" = 0.65), guide = "none") +
  scale_x_log10(breaks = c(100, 250, 500, 1000, 2500)) +
  scale_color_manual(values = pal) +
  facet_wrap(~ Group_Name, nrow = 1) +
  coord_cartesian(ylim = c(0.15, 1.0)) +
  labs(
    title = "Empirical Selection Consistency & Noise Robustness Across Dynamical Regimes",
    subtitle = "Top: F1 vs. Sample Size N (marginalized over noise) | Bottom: F1 vs. Noise Level sigma (marginalized over N)",
    x = "Sample Size N (log scale)",
    y = "Structure Recovery F1",
    color = "Method"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "#f0f4f8"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

p_noise <- ggplot(df_noise, aes(x = Noise, y = Mean_F1, color = Method_Name, group = Method_Name)) +
  geom_line(aes(linewidth = Method_Name == "BGNLM (Tuned)"), alpha = 0.92) +
  geom_point(size = 1.8) +
  scale_linewidth_manual(values = c("TRUE" = 1.35, "FALSE" = 0.65), guide = "none") +
  scale_x_log10(breaks = c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0),
                labels = c("0.001", "0.01", "0.05", "0.1", "0.25", "0.5", "1.0")) +
  scale_color_manual(values = pal) +
  facet_wrap(~ Group_Name, nrow = 1) +
  coord_cartesian(ylim = c(0.15, 1.0)) +
  labs(
    x = "Noise Intensity sigma (log scale)",
    y = "Structure Recovery F1",
    color = "Method"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 8.5),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

# Stack p_n and p_noise cleanly using gridExtra / gtable / patchwork fallback
if (requireNamespace("gridExtra", quietly = TRUE)) {
  g <- gridExtra::arrangeGrob(p_n, p_noise, ncol = 1, heights = c(1.0, 1.15))
  ggplot2::ggsave("results/figures/fig_consistency_scaling.pdf", g, width = 11, height = 5.4, bg = "white")
  ggplot2::ggsave("results/figures/fig_consistency_scaling.png", g, width = 11, height = 5.4, dpi = 200, bg = "white")
  ggplot2::ggsave("paper_icml/figures/fig_consistency_scaling.pdf", g, width = 11, height = 5.4, bg = "white")
} else {
  # Pure grid viewport composition
  save_stacked <- function(path_out, dev_fn) {
    dev_fn(path_out)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1, heights = grid::unit(c(0.47, 0.53), "npc"))))
    print(p_n, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_noise, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grDevices::dev.off()
  }
  save_stacked("results/figures/fig_consistency_scaling.pdf", function(p) grDevices::pdf(p, width = 11, height = 5.4))
  save_stacked("paper_icml/figures/fig_consistency_scaling.pdf", function(p) grDevices::pdf(p, width = 11, height = 5.4))
  save_stacked("results/figures/fig_consistency_scaling.png", function(p) grDevices::png(p, width = 11, height = 5.4, units = "in", res = 200))
}

# ------------------------------------------------------------------------------
# Figure B: fig_uq_stability_pareto.pdf (Synthetic Stability vs F1 & Real Parsimony vs F1)
# ------------------------------------------------------------------------------
all_short <- c(
  "BGNLM_Full_Tuned"   = "BGNLM (Tuned)",
  "BGLM_Oracle_Basis"  = "BGLM (Oracle)",
  "BGNLM_Hybrid_Poly"  = "BGNLM-Hybrid",
  "BGLM_Poly_Basis"    = "BGLM (Poly)",
  "SINDy_STLSQ_Oracle" = "STLSQ (Oracle)",
  "BGLM_Overcomplete"  = "BGLM (Overcomplete)",
  "SINDy_STLSQ_Poly"   = "STLSQ (Poly)",
  "BGNLM_IWSM_Untuned" = "BGNLM (Untuned)",
  "UQ_SINDy_Oracle"    = "Bag-STLSQ (Oracle)",
  "UQ_SINDy_Poly"      = "Bag-STLSQ (Poly)",
  "SINDy_PI_Rational"  = "SINDy-RatPoly2",
  "SR_Operon"          = "SR Pool (MDL)",
  "SR_PySR"            = "SR Pool (Pareto)"
)

fam_map <- c(
  "BGNLM_Full_Tuned"   = "BGNLM (Basis-Free)",
  "BGNLM_Hybrid_Poly"  = "BGNLM (Basis-Free)",
  "BGNLM_IWSM_Untuned" = "BGNLM (Basis-Free)",
  "BGLM_Oracle_Basis"  = "Predefined BGLM",
  "BGLM_Poly_Basis"    = "Predefined BGLM",
  "BGLM_Overcomplete"  = "Predefined BGLM",
  "SINDy_STLSQ_Oracle" = "SINDy / Bag-STLSQ",
  "SINDy_STLSQ_Poly"   = "SINDy / Bag-STLSQ",
  "UQ_SINDy_Oracle"    = "SINDy / Bag-STLSQ",
  "UQ_SINDy_Poly"      = "SINDy / Bag-STLSQ",
  "SINDy_PI_Rational"  = "SINDy / Bag-STLSQ",
  "SR_Operon"          = "Static-Pool SR",
  "SR_PySR"            = "Static-Pool SR"
)

fam_pal <- c(
  "BGNLM (Basis-Free)"  = "#0072B2",
  "Predefined BGLM"     = "#009E73",
  "SINDy / Bag-STLSQ"   = "#D55E00",
  "Static-Pool SR"      = "#9467BD"
)

synth_m <- synth_stab[, .(
  Mean_F1 = mean(Mean_F1, na.rm = TRUE),
  Mean_Nogueira = mean(Nogueira_Phi, na.rm = TRUE),
  Mean_FDR = mean(Mean_FDR, na.rm = TRUE)
), by = Method]
synth_m[, Label := all_short[Method]]
synth_m[, Family := fam_map[Method]]

p_stab <- ggplot(synth_m, aes(x = Mean_Nogueira, y = Mean_F1, color = Family)) +
  geom_point(aes(size = 1 - Mean_FDR), alpha = 0.88) +
  geom_text(aes(label = Label), size = 2.6, vjust = -0.85, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = fam_pal) +
  scale_size_continuous(range = c(2.5, 6.0), guide = "none") +
  xlim(0.42, 0.96) + ylim(0.38, 0.94) +
  labs(
    title = "(A) Synthetic Suite: Stability vs. Structure F1",
    x = expression("Nogueira Selection Stability " * hat(Phi) * " (higher is better)"),
    y = "Mean Structure Recovery F1",
    color = "Method Family"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 10)
  )

real_m <- real_stab[, .(
  Mean_F1 = mean(Mean_F1, na.rm = TRUE),
  Mean_OC = mean(Mean_OC, na.rm = TRUE),
  Mean_R2_traj = mean(Mean_R2_traj, na.rm = TRUE)
), by = Method]
real_m[, Label := all_short[Method]]
real_m[, Family := fam_map[Method]]

p_real <- ggplot(real_m, aes(x = Mean_OC, y = Mean_F1, color = Family)) +
  geom_point(aes(size = pmax(0.2, Mean_R2_traj)), alpha = 0.88) +
  geom_text(aes(label = Label), size = 2.6, vjust = -0.85, hjust = 0.5, show.legend = FALSE) +
  scale_color_manual(values = fam_pal) +
  scale_size_continuous(range = c(2.5, 6.0), name = expression("Rollout " * R[traj]^2)) +
  xlim(0, 62) + ylim(0.22, 0.78) +
  labs(
    title = "(B) Quasi-Empirical Suite: Parsimony (OC) vs. F1",
    x = "Mean Operation Count OC (lower is simpler)",
    y = "Mean Structure Recovery F1",
    color = "Method Family"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 10)
  ) +
  guides(size = "none")

if (requireNamespace("gridExtra", quietly = TRUE)) {
  g2 <- gridExtra::arrangeGrob(p_stab, p_real, ncol = 2)
  ggplot2::ggsave("results/figures/fig_uq_stability_pareto.pdf", g2, width = 11.2, height = 4.6, bg = "white")
  ggplot2::ggsave("results/figures/fig_uq_stability_pareto.png", g2, width = 11.2, height = 4.6, dpi = 200, bg = "white")
  ggplot2::ggsave("paper_icml/figures/fig_uq_stability_pareto.pdf", g2, width = 11.2, height = 4.6, bg = "white")
} else {
  save_side <- function(path_out, dev_fn) {
    dev_fn(path_out)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
    print(p_stab, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(p_real, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    grDevices::dev.off()
  }
  save_side("results/figures/fig_uq_stability_pareto.pdf", function(p) grDevices::pdf(p, width = 11.2, height = 4.6))
  save_side("paper_icml/figures/fig_uq_stability_pareto.pdf", function(p) grDevices::pdf(p, width = 11.2, height = 4.6))
  save_side("results/figures/fig_uq_stability_pareto.png", function(p) grDevices::png(p, width = 11.2, height = 4.6, units = "in", res = 200))
}

message("[OK] Generated fig_consistency_scaling and fig_uq_stability_pareto in results/figures/ and paper_icml/figures/.")
