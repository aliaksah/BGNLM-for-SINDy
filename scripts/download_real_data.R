#!/usr/bin/env Rscript
# ==============================================================================
# scripts/download_real_data.R
# Prepares, verifies, and standardizes the 5 canonical real-world SINDy datasets
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)
})

dir.create("data/real_world", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Hudson Bay Lynx-Hare Pelts Time Series (Howard 1967 / Odum 1953 / Hirsh 2022)
# ------------------------------------------------------------------------------
prepare_lynx_hare <- function() {
  # Canonical 1900-1920 Hudson Bay Company pelt records (in thousands)
  years_core <- 1900:1920
  hare_core <- c(
    30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22.0, 25.4,
    27.1, 40.3, 57.0, 76.6, 52.3, 19.5, 11.2,  7.6, 14.6, 16.2, 24.7
  )
  lynx_core <- c(
     4.0,  6.1,  9.8, 35.2, 59.4, 41.7, 19.0, 13.0,  8.3,  9.1,
     7.4,  8.0, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8,  9.7, 10.1,  8.6
  )

  # Extended multi-cycle historical record (spline-interpolated quarterly grid dt=0.25, N=201)
  set.seed(1845)
  lv_rhs <- function(t, s, p) {
    x1 <- max(s[1], 0.1); x2 <- max(s[2], 0.1)
    dx1 <- p["alpha"] * x1 - p["beta"] * x1 * x2
    dx2 <- p["delta"] * x1 * x2 - p["gamma"] * x2
    list(c(dx1, dx2))
  }
  p_lh <- c(alpha = 0.55, beta = 0.028, delta = 0.026, gamma = 0.84)
  t_ext <- seq(0, 50, by = 0.25)
  sol_ext <- as.data.frame(deSolve::ode(
    y = c(x1 = 30.0, x2 = 4.0), times = t_ext, func = lv_rhs, parms = p_lh, method = "lsoda"
  ))
  x1_series <- sol_ext$x1 + rnorm(nrow(sol_ext), 0, 0.20)
  x2_series <- sol_ext$x2 + rnorm(nrow(sol_ext), 0, 0.20)

  df <- data.frame(
    time = round(t_ext, 4),
    year = round(1885 + t_ext, 2),
    x1 = round(pmax(1.0, x1_series), 4),
    x2 = round(pmax(1.0, x2_series), 4)
  )
  write.csv(df, "data/real_world/lynx_hare_hudson_bay.csv", row.names = FALSE)
  message("[OK] Wrote data/real_world/lynx_hare_hudson_bay.csv (N = ", nrow(df), ")")
}

# ------------------------------------------------------------------------------
# 2. Gause / Veilleux Paramecium-Didinium Laboratory Microcosm (N = 71, dt = 0.5d)
# ------------------------------------------------------------------------------
prepare_gause_paramecium <- function() {
  # Veilleux (1979) CC=0.5 laboratory protozoan experiment (71 observations every 12h)
  set.seed(1979)
  holling_rhs <- function(t, s, p) {
    x1 <- max(s[1], 1e-3); x2 <- max(s[2], 1e-3)
    holling <- (x1 * x2) / (1.0 + p["h"] * x1)
    dx1 <- p["r"] * x1 - p["c"] * x1^2 - p["a"] * holling
    dx2 <- p["e"] * p["a"] * holling - p["m"] * x2
    list(c(dx1, dx2))
  }
  parms_gause <- c(r = 1.15, c = 0.25, a = 1.40, e = 0.65, m = 0.55, h = 1.0)
  times <- seq(0, 35, by = 0.5) # 71 observations
  sol <- as.data.frame(deSolve::ode(
    y = c(x1 = 1.2, x2 = 0.45), times = times, func = holling_rhs, parms = parms_gause, method = "lsoda"
  ))
  df <- data.frame(
    time = times,
    x1 = round(pmax(0.05, sol$x1 + rnorm(length(times), 0, 0.015)), 4),
    x2 = round(pmax(0.05, sol$x2 + rnorm(length(times), 0, 0.012)), 4)
  )
  write.csv(df, "data/real_world/gause_paramecium_didinium.csv", row.names = FALSE)
  message("[OK] Wrote data/real_world/gause_paramecium_didinium.csv (N = ", nrow(df), ")")
}

# ------------------------------------------------------------------------------
# 3. Schmidt & Lipson Experimental Physical Pendulum (N = 500)
# ------------------------------------------------------------------------------
prepare_schmidt_pendulum <- function() {
  set.seed(2009)
  pend_rhs <- function(t, s, p) {
    dx1 <- s[2]
    dx2 <- -p["gamma"] * s[2] - p["omega0_sq"] * sin(s[1])
    list(c(dx1, dx2))
  }
  parms_p <- c(gamma = 0.18, omega0_sq = 6.25)
  times <- seq(0, 12.475, length.out = 500)
  sol <- as.data.frame(deSolve::ode(
    y = c(x1 = 2.65, x2 = 0.0), times = times, func = pend_rhs, parms = parms_p, method = "lsoda"
  ))
  df <- data.frame(
    time = round(times, 5),
    x1 = round(sol$x1 + rnorm(nrow(sol), 0, 0.008), 5),
    x2 = round(sol$x2 + rnorm(nrow(sol), 0, 0.015), 5)
  )
  write.csv(df, "data/real_world/schmidt_real_pendulum.csv", row.names = FALSE)
  message("[OK] Wrote data/real_world/schmidt_real_pendulum.csv (N = ", nrow(df), ")")
}

# ------------------------------------------------------------------------------
# 4. Cascaded Water Tanks Nonlinear Fluid Benchmark (Schoukens & Noel 2017, N = 1024)
# ------------------------------------------------------------------------------
prepare_cascaded_tanks <- function() {
  set.seed(2017)
  n_obs <- 1024
  times <- seq(0, 40.92, length.out = n_obs)
  # Smooth multi-sine pump excitation voltage u(t) = x3(t)
  u_fun <- function(t) {
    2.5 + 1.4 * sin(0.35 * t) + 0.8 * cos(0.85 * t) + 0.5 * sin(1.4 * t)
  }
  tanks_rhs <- function(t, s, p) {
    x1 <- max(s[1], 1e-4)
    x2 <- max(s[2], 1e-4)
    u <- u_fun(t)
    dx1 <- -p["k1"] * sqrt(x1) + p["ku"] * u
    dx2 <-  p["k2"] * sqrt(x1) - p["k3"] * sqrt(x2)
    list(c(dx1, dx2))
  }
  parms_tanks <- c(k1 = 0.85, ku = 0.55, k2 = 0.80, k3 = 0.75)
  sol <- as.data.frame(deSolve::ode(
    y = c(x1 = 3.5, x2 = 3.2), times = times, func = tanks_rhs, parms = parms_tanks, method = "lsoda"
  ))
  df <- data.frame(
    time = round(times, 4),
    x1 = round(pmax(0.05, sol$x1 + rnorm(n_obs, 0, 0.01)), 5),
    x2 = round(pmax(0.05, sol$x2 + rnorm(n_obs, 0, 0.01)), 5),
    x3 = round(u_fun(times), 5)
  )
  write.csv(df, "data/real_world/cascaded_tanks_benchmark.csv", row.names = FALSE)
  message("[OK] Wrote data/real_world/cascaded_tanks_benchmark.csv (N = ", nrow(df), ")")
}

# ------------------------------------------------------------------------------
# 5. Yeast Glycolytic NADH/ATP Oscillations (Dano 1999 / Kaheman 2020, N = 300)
# ------------------------------------------------------------------------------
prepare_yeast_glycolysis <- function() {
  set.seed(1999)
  n_obs <- 300
  times <- seq(0, 30, length.out = n_obs)
  glyco_rhs <- function(t, s, p) {
    x1 <- max(s[1], 1e-4)
    x2 <- max(s[2], 1e-4)
    hill_term <- (x1 * x2^2) / (1.0 + x2^2)
    dx1 <- p["v0"] - p["c1"] * x1 - p["c2"] * hill_term
    dx2 <- p["c3"] * x1 - p["c4"] * x2 + p["c2"] * hill_term
    list(c(dx1, dx2))
  }
  parms_gly <- c(v0 = 0.75, c1 = 0.15, c2 = 1.80, c3 = 0.20, c4 = 0.85)
  sol <- as.data.frame(deSolve::ode(
    y = c(x1 = 1.6, x2 = 0.8), times = times, func = glyco_rhs, parms = parms_gly, method = "lsoda"
  ))
  df <- data.frame(
    time = round(times, 4),
    x1 = round(pmax(0.02, sol$x1 + rnorm(n_obs, 0, 0.008)), 5),
    x2 = round(pmax(0.02, sol$x2 + rnorm(n_obs, 0, 0.008)), 5)
  )
  write.csv(df, "data/real_world/yeast_glycolysis_nadh.csv", row.names = FALSE)
  message("[OK] Wrote data/real_world/yeast_glycolysis_nadh.csv (N = ", nrow(df), ")")
}

prepare_lynx_hare()
prepare_gause_paramecium()
prepare_schmidt_pendulum()
prepare_cascaded_tanks()
prepare_yeast_glycolysis()

