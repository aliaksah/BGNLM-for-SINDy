# ==============================================================================
# R/01_ode_simulators.R
# 12 Canonical SINDy Benchmark Dynamical Systems, Noise Regimes, & Differentiation
# ==============================================================================

suppressPackageStartupMessages({
  library(deSolve)
  library(signal)
})

#' Return complete Registry of the 12 Benchmark Dynamical Systems
get_systems_registry <- function() {
  list(
    sys01_linear_quadratic_3d = list(
      sys_id = 1L,
      id = "sys01_linear_quadratic_3d",
      name = "Linear-Quadratic Coupled 3D",
      group = "Group_A_Polynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = 2.0, x2 = 0.0, x3 = 1.0),
      parameters = c(a = -1.0, b = 20.0, c = -20.0, d = -1.0, e = -3.0),
      dt = 0.0001,
      T_max = 20.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- parms["a"] * x1 + parms["b"] * x1 * x2
        dx2 <- parms["c"] * x1 + parms["d"] * x1 * x2
        dx3 <- parms["e"] * x3
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = parms["a"] * x1 + parms["b"] * x1 * x2,
          dx2 = parms["c"] * x1 + parms["d"] * x1 * x2,
          dx3 = parms["e"] * x3
        )
      },
      ground_truth = list(
        c("x1", "x1*x2"),
        c("x1", "x1*x2"),
        c("x3")
      ),
      true_coefs = list(
        c("x1" = -1.0, "x1*x2" = 20.0),
        c("x1" = -20.0, "x1*x2" = -1.0),
        c("x3" = -3.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys02_lorenz63 = list(
      sys_id = 2L,
      id = "sys02_lorenz63",
      name = "Lorenz-63 System",
      group = "Group_A_Polynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = -0.5, x2 = -2.0, x3 = 3.0),
      parameters = c(sigma = 10.0, rho = 28.0, beta = 8.0 / 3.0),
      dt = 0.0001,
      T_max = 50.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- parms["sigma"] * (x2 - x1)
        dx2 <- x1 * (parms["rho"] - x3) - x2
        dx3 <- x1 * x2 - parms["beta"] * x3
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = parms["sigma"] * (x2 - x1),
          dx2 = x1 * (parms["rho"] - x3) - x2,
          dx3 = x1 * x2 - parms["beta"] * x3
        )
      },
      ground_truth = list(
        c("x1", "x2"),
        c("x1", "x2", "x1*x3"),
        c("x3", "x1*x2")
      ),
      true_coefs = list(
        c("x1" = -10.0, "x2" = 10.0),
        c("x1" = 28.0, "x2" = -1.0, "x1*x3" = -1.0),
        c("x3" = -8.0 / 3.0, "x1*x2" = 1.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys03_rossler = list(
      sys_id = 3L,
      id = "sys03_rossler",
      name = "Rossler Attractor",
      group = "Group_A_Polynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = 0.5, x2 = -1.0, x3 = 2.0),
      parameters = c(a = 0.2, b = 0.2, c = 5.7),
      dt = 0.0001,
      T_max = 50.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- -x2 - x3
        dx2 <- x1 + parms["a"] * x2
        dx3 <- parms["b"] + x3 * (x1 - parms["c"])
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = -x2 - x3,
          dx2 = x1 + parms["a"] * x2,
          dx3 = parms["b"] + x3 * (x1 - parms["c"])
        )
      },
      ground_truth = list(
        c("x2", "x3"),
        c("x1", "x2"),
        c("x3", "x1*x3")
      ),
      true_coefs = list(
        c("x2" = -1.0, "x3" = -1.0),
        c("x1" = 1.0, "x2" = 0.2),
        c("x3" = -5.7, "x1*x3" = 1.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys04_lotka_volterra = list(
      sys_id = 4L,
      id = "sys04_lotka_volterra",
      name = "Lotka-Volterra Predator-Prey",
      group = "Group_A_Polynomial",
      dim = 2L,
      state_names = c("x1", "x2"),
      initial_conditions = c(x1 = 1.0, x2 = 1.0),
      parameters = c(alpha = 1.5, beta = 1.0, delta = 1.0, gamma = 3.0),
      dt = 0.0001,
      T_max = 30.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]
        dx1 <- parms["alpha"] * x1 - parms["beta"] * x1 * x2
        dx2 <- parms["delta"] * x1 * x2 - parms["gamma"] * x2
        list(c(dx1, dx2))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]
        cbind(
          dx1 = parms["alpha"] * x1 - parms["beta"] * x1 * x2,
          dx2 = parms["delta"] * x1 * x2 - parms["gamma"] * x2
        )
      },
      ground_truth = list(
        c("x1", "x1*x2"),
        c("x2", "x1*x2")
      ),
      true_coefs = list(
        c("x1" = 1.5, "x1*x2" = -1.0),
        c("x2" = -3.0, "x1*x2" = 1.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys05_van_der_pol = list(
      sys_id = 5L,
      id = "sys05_van_der_pol",
      name = "Van der Pol Oscillator",
      group = "Group_A_Polynomial",
      dim = 2L,
      state_names = c("x1", "x2"),
      initial_conditions = c(x1 = 2.0, x2 = 0.0),
      parameters = c(mu = 2.0),
      dt = 0.0001,
      T_max = 30.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]
        dx1 <- x2
        dx2 <- parms["mu"] * (1.0 - x1^2) * x2 - x1
        list(c(dx1, dx2))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]
        cbind(
          dx1 = x2,
          dx2 = parms["mu"] * (1.0 - x1^2) * x2 - x1
        )
      },
      ground_truth = list(
        c("x2"),
        c("x1", "x2", "x1^2*x2")
      ),
      true_coefs = list(
        c("x2" = 1.0),
        c("x1" = -1.0, "x2" = 2.0, "x1^2*x2" = -2.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys06_duffing = list(
      sys_id = 6L,
      id = "sys06_duffing",
      name = "Duffing Oscillator",
      group = "Group_A_Polynomial",
      dim = 2L,
      state_names = c("x1", "x2"),
      initial_conditions = c(x1 = 0.5, x2 = 1.0),
      parameters = c(delta = 0.1, alpha = -1.0, beta = 1.0),
      dt = 0.0001,
      T_max = 30.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]
        dx1 <- x2
        dx2 <- -parms["delta"] * x2 - parms["alpha"] * x1 - parms["beta"] * x1^3
        list(c(dx1, dx2))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]
        cbind(
          dx1 = x2,
          dx2 = -parms["delta"] * x2 - parms["alpha"] * x1 - parms["beta"] * x1^3
        )
      },
      ground_truth = list(
        c("x2"),
        c("x1", "x2", "x1^3")
      ),
      true_coefs = list(
        c("x2" = 1.0),
        c("x1" = 1.0, "x2" = -0.1, "x1^3" = -1.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys07_brusselator = list(
      sys_id = 7L,
      id = "sys07_brusselator",
      name = "Brusselator System",
      group = "Group_A_Polynomial",
      dim = 2L,
      state_names = c("x1", "x2"),
      initial_conditions = c(x1 = 1.5, x2 = 3.0),
      parameters = c(a = 1.0, b = 3.0),
      dt = 0.0001,
      T_max = 30.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]
        dx1 <- parms["a"] - (parms["b"] + 1.0) * x1 + (x1^2) * x2
        dx2 <- parms["b"] * x1 - (x1^2) * x2
        list(c(dx1, dx2))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]
        cbind(
          dx1 = parms["a"] - (parms["b"] + 1.0) * x1 + (x1^2) * x2,
          dx2 = parms["b"] * x1 - (x1^2) * x2
        )
      },
      ground_truth = list(
        c("x1", "x1^2*x2"),
        c("x1", "x1^2*x2")
      ),
      true_coefs = list(
        c("x1" = -4.0, "x1^2*x2" = 1.0),
        c("x1" = 3.0, "x1^2*x2" = -1.0)
      ),
      oracle_extra_basis = character(0)
    ),

    sys08_rossler_lorenz_hybrid = list(
      sys_id = 8L,
      id = "sys08_rossler_lorenz_hybrid",
      name = "Rossler-Lorenz Hybrid 3D",
      group = "Group_B_NonPolynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = 0.5, x2 = -1.0, x3 = 2.0),
      parameters = c(a = 0.2, b = 0.1, c = 0.4, d = 5.7),
      dt = 0.0001,
      T_max = 50.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- -x2 - x3 + parms["a"] * sin(x1)
        dx2 <- x1 + parms["b"] * x2
        dx3 <- parms["c"] + x3 * (x1 - parms["d"])
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = -x2 - x3 + parms["a"] * sin(x1),
          dx2 = x1 + parms["b"] * x2,
          dx3 = parms["c"] + x3 * (x1 - parms["d"])
        )
      },
      ground_truth = list(
        c("x2", "x3", "sin(x1)"),
        c("x1", "x2"),
        c("x3", "x1*x3")
      ),
      true_coefs = list(
        c("x2" = -1.0, "x3" = -1.0, "sin(x1)" = 0.2),
        c("x1" = 1.0, "x2" = 0.1),
        c("x3" = -5.7, "x1*x3" = 1.0)
      ),
      oracle_extra_basis = c("sin(x1)", "cos(x1)")
    ),

    sys09_damped_pendulum = list(
      sys_id = 9L,
      id = "sys09_damped_pendulum",
      name = "Damped Nonlinear Pendulum",
      group = "Group_B_NonPolynomial",
      dim = 2L,
      state_names = c("x1", "x2"),
      initial_conditions = c(x1 = pi - 0.1, x2 = 0.0),
      parameters = c(gamma = 0.25, g_over_L = 5.0),
      dt = 0.0001,
      T_max = 25.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]
        dx1 <- x2
        dx2 <- -parms["gamma"] * x2 - parms["g_over_L"] * sin(x1)
        list(c(dx1, dx2))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]
        cbind(
          dx1 = x2,
          dx2 = -parms["gamma"] * x2 - parms["g_over_L"] * sin(x1)
        )
      },
      ground_truth = list(
        c("x2"),
        c("x2", "sin(x1)")
      ),
      true_coefs = list(
        c("x2" = 1.0),
        c("x2" = -0.25, "sin(x1)" = -5.0)
      ),
      oracle_extra_basis = c("sin(x1)", "cos(x1)")
    ),

    sys10_thomas_attractor = list(
      sys_id = 10L,
      id = "sys10_thomas_attractor",
      name = "Thomas Cyclically Symmetric Attractor",
      group = "Group_B_NonPolynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = 0.1, x2 = 0.0, x3 = 0.0),
      parameters = c(b = 0.18),
      dt = 0.0001,
      T_max = 50.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- sin(x2) - parms["b"] * x1
        dx2 <- sin(x3) - parms["b"] * x2
        dx3 <- sin(x1) - parms["b"] * x3
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = sin(x2) - parms["b"] * x1,
          dx2 = sin(x3) - parms["b"] * x2,
          dx3 = sin(x1) - parms["b"] * x3
        )
      },
      ground_truth = list(
        c("x1", "sin(x2)"),
        c("x2", "sin(x3)"),
        c("x3", "sin(x1)")
      ),
      true_coefs = list(
        c("x1" = -0.18, "sin(x2)" = 1.0),
        c("x2" = -0.18, "sin(x3)" = 1.0),
        c("x3" = -0.18, "sin(x1)" = 1.0)
      ),
      oracle_extra_basis = c("sin(x1)", "sin(x2)", "sin(x3)", "cos(x1)", "cos(x2)", "cos(x3)")
    ),

    sys11_goodwin_mm = list(
      sys_id = 11L,
      id = "sys11_goodwin_mm",
      name = "Michaelis-Menten Goodwin Oscillator",
      group = "Group_B_NonPolynomial",
      dim = 3L,
      state_names = c("x1", "x2", "x3"),
      initial_conditions = c(x1 = 0.5, x2 = 0.5, x3 = 0.5),
      parameters = c(p = 2.0, v0 = 1.0, k1 = 0.5, k2 = 0.5, vm = 1.0, Km = 1.0),
      dt = 0.0001,
      T_max = 40.0,
      rhs = function(t, state, parms) {
        x1 <- state[1]; x2 <- state[2]; x3 <- state[3]
        dx1 <- parms["v0"] / (1.0 + x3^2) - parms["k1"] * x1
        dx2 <- parms["k1"] * x1 - parms["k2"] * x2
        dx3 <- parms["k2"] * x2 - (parms["vm"] * x3) / (parms["Km"] + abs(x3))
        list(c(dx1, dx2, dx3))
      },
      analytical_deriv = function(X, parms) {
        x1 <- X[, 1]; x2 <- X[, 2]; x3 <- X[, 3]
        cbind(
          dx1 = parms["v0"] / (1.0 + x3^2) - parms["k1"] * x1,
          dx2 = parms["k1"] * x1 - parms["k2"] * x2,
          dx3 = parms["k2"] * x2 - (parms["vm"] * x3) / (parms["Km"] + abs(x3))
        )
      },
      ground_truth = list(
        c("x1", "1/(1+x3^2)"),
        c("x1", "x2"),
        c("x2", "x3/(1+x3)")
      ),
      true_coefs = list(
        c("x1" = -0.5, "1/(1+x3^2)" = 1.0),
        c("x1" = 0.5, "x2" = -0.5),
        c("x2" = 0.5, "x3/(1+x3)" = -1.0)
      ),
      oracle_extra_basis = c("1/(1+x1^2)", "1/(1+x2^2)", "1/(1+x3^2)", "x1/(1+x1)", "x2/(1+x2)", "x3/(1+x3)")
    ),

    sys12_lorenz96_5d = list(
      sys_id = 12L,
      id = "sys12_lorenz96_5d",
      name = "Lorenz-96 (5D Coupled)",
      group = "Group_C_HighDim",
      dim = 5L,
      state_names = c("x1", "x2", "x3", "x4", "x5"),
      initial_conditions = c(x1 = 8.0, x2 = 8.0, x3 = 8.0, x4 = 8.0, x5 = 8.01),
      parameters = c(F = 8.0),
      dt = 0.0001,
      T_max = 30.0,
      rhs = function(t, state, parms) {
        x <- state
        F_val <- parms["F"]
        dx1 <- (x[2] - x[4]) * x[5] - x[1] + F_val
        dx2 <- (x[3] - x[5]) * x[1] - x[2] + F_val
        dx3 <- (x[4] - x[1]) * x[2] - x[3] + F_val
        dx4 <- (x[5] - x[2]) * x[3] - x[4] + F_val
        dx5 <- (x[1] - x[3]) * x[4] - x[5] + F_val
        list(c(dx1, dx2, dx3, dx4, dx5))
      },
      analytical_deriv = function(X, parms) {
        F_val <- parms["F"]
        cbind(
          dx1 = (X[, 2] - X[, 4]) * X[, 5] - X[, 1] + F_val,
          dx2 = (X[, 3] - X[, 5]) * X[, 1] - X[, 2] + F_val,
          dx3 = (X[, 4] - X[, 1]) * X[, 2] - X[, 3] + F_val,
          dx4 = (X[, 5] - X[, 2]) * X[, 3] - X[, 4] + F_val,
          dx5 = (X[, 1] - X[, 3]) * X[, 4] - X[, 5] + F_val
        )
      },
      ground_truth = list(
        c("x1", "x2*x5", "x4*x5"),
        c("x2", "x1*x3", "x1*x5"),
        c("x3", "x1*x2", "x2*x4"),
        c("x4", "x2*x3", "x3*x5"),
        c("x5", "x1*x4", "x3*x4")
      ),
      true_coefs = list(
        c("x1" = -1.0, "x2*x5" = 1.0, "x4*x5" = -1.0),
        c("x2" = -1.0, "x1*x3" = 1.0, "x1*x5" = -1.0),
        c("x3" = -1.0, "x1*x2" = -1.0, "x2*x4" = 1.0),
        c("x4" = -1.0, "x2*x3" = -1.0, "x3*x5" = 1.0),
        c("x5" = -1.0, "x1*x4" = 1.0, "x3*x4" = -1.0)
      ),
      oracle_extra_basis = character(0)
    )
  )
}

#' Finite difference derivative helper (compatible with original IWSM scripts)
finite_difference_matrix <- function(X, dt) {
  dX <- apply(X, 2, diff) / dt
  dX
}

#' Savitzky-Golay + Smoothing Spline numerical derivative estimator
estimate_derivatives_sgolay_spline <- function(X, times, sg_p = 3, sg_n = 7) {
  dt <- median(diff(times))
  n_obs <- nrow(X)
  m <- ncol(X)
  dX <- matrix(0, nrow = n_obs, ncol = m)
  colnames(dX) <- paste0("d", colnames(X))

  for (j in seq_len(m)) {
    y_col <- X[, j]
    sp_fit <- tryCatch(
      stats::smooth.spline(x = times, y = y_col, spar = 0.30),
      error = function(e) NULL
    )
    if (!is.null(sp_fit)) {
      d_sp <- stats::predict(sp_fit, x = times, deriv = 1)$y
      if (n_obs >= sg_n && sg_n %% 2 == 1) {
        d_sg <- tryCatch(
          signal::sgolayfilt(y_col, p = sg_p, n = sg_n, m = 1) / dt,
          error = function(e) d_sp
        )
        dX[, j] <- 0.5 * (d_sp + d_sg)
      } else {
        dX[, j] <- d_sp
      }
    } else {
      # Central difference fallback
      d_col <- numeric(n_obs)
      d_col[1] <- (y_col[2] - y_col[1]) / dt
      d_col[n_obs] <- (y_col[n_obs] - y_col[n_obs - 1]) / dt
      if (n_obs > 2) {
        d_col[2:(n_obs - 1)] <- (y_col[3:n_obs] - y_col[1:(n_obs - 2)]) / (2 * dt)
      }
      dX[, j] <- d_col
    }
  }
  dX
}

#' Total Variation Regularized Differentiation (Chartrand 2011 TVRegDiff in R)
estimate_derivatives_tvregdiff <- function(X, times, alpha = 1e-2, n_iter = 15) {
  dt <- median(diff(times))
  n_obs <- nrow(X)
  m <- ncol(X)
  dX <- matrix(0, nrow = n_obs, ncol = m)
  colnames(dX) <- paste0("d", colnames(X))

  for (j in seq_len(m)) {
    u <- X[, j]
    # Initialize with central difference
    du <- numeric(n_obs)
    du[1] <- (u[2] - u[1]) / dt
    du[n_obs] <- (u[n_obs] - u[n_obs - 1]) / dt
    if (n_obs > 2) {
      du[2:(n_obs - 1)] <- (u[3:n_obs] - u[1:(n_obs - 2)]) / (2 * dt)
    }
    # Iterative TV-penalized tridiagonal smoothing of derivative
    eps_tv <- 1e-4 * max(sd(du), 1e-6)
    for (iter in seq_len(n_iter)) {
      d2 <- diff(du)
      w <- 1.0 / sqrt(d2^2 + eps_tv^2)
      # Build penalized normal equations: (I + alpha * D^T W D) du_new = du_init
      diag_main <- rep(1.0, n_obs)
      diag_main[1:(n_obs - 1)] <- diag_main[1:(n_obs - 1)] + alpha * w
      diag_main[2:n_obs] <- diag_main[2:n_obs] + alpha * w
      off_diag <- -alpha * w
      M <- diag(diag_main)
      M[cbind(1:(n_obs - 1), 2:n_obs)] <- off_diag
      M[cbind(2:n_obs, 1:(n_obs - 1))] <- off_diag
      du <- tryCatch(solve(M, du), error = function(e) du)
    }
    dX[, j] <- du
  }
  dX
}

# In-memory cache environment for base noiseless trajectories
.traj_cache <- new.env(parent = emptyenv())

#' Simulate or retrieve cached fine-grid trajectory for a synthetic ODE system
simulate_base_trajectory <- function(system_id, dt = NULL, T_max = NULL, use_cache = TRUE) {
  reg <- get_systems_registry()
  if (!system_id %in% names(reg)) {
    stop("Unknown system_id: ", system_id)
  }
  sys <- reg[[system_id]]
  dt_use <- if (is.null(dt)) sys$dt else dt
  T_use <- if (is.null(T_max)) sys$T_max else T_max
  cache_key <- paste0(system_id, "_dt", dt_use, "_T", T_use)

  if (use_cache && exists(cache_key, envir = .traj_cache)) {
    return(get(cache_key, envir = .traj_cache))
  }

  times <- seq(0, T_use, by = dt_use)
  sol <- deSolve::ode(
    y = sys$initial_conditions,
    times = times,
    func = sys$rhs,
    parms = sys$parameters,
    method = "lsoda",
    rtol = 1e-9,
    atol = 1e-9
  )
  sol_mat <- as.matrix(sol)
  X_full <- sol_mat[, -1, drop = FALSE]
  colnames(X_full) <- sys$state_names

  # Finite difference derivative (aligned with X_full[-1, ])
  dX_fd <- finite_difference_matrix(X_full, dt_use)
  # Analytical derivative at X_full[-1, ]
  X_aligned <- X_full[-1, , drop = FALSE]
  t_aligned <- sol_mat[-1, 1]
  dX_analytical <- sys$analytical_deriv(X_aligned, sys$parameters)
  colnames(dX_analytical) <- paste0("d", sys$state_names)
  colnames(dX_fd) <- paste0("d", sys$state_names)

  out <- list(
    system = sys,
    times = t_aligned,
    X = X_aligned,
    dX_analytical = dX_analytical,
    dX_fd = dX_fd,
    dt = dt_use,
    T_max = T_use
  )
  if (use_cache) {
    assign(cache_key, out, envir = .traj_cache)
  }
  out
}

#' Generate a complete synthetic benchmark task split (Train, In-Sample Test, OOD Test)
generate_synthetic_task_data <- function(
  system_id,
  N_train = 500,
  noise_eta = 0.05,
  noise_regime = c("derivative_snr", "state_noise", "iwsm_absolute"),
  seed = 10101,
  N_test_in = 1000,
  N_test_ood = 1000,
  dt = NULL,
  T_max = NULL
) {
  noise_regime <- match.arg(noise_regime)
  base <- simulate_base_trajectory(system_id, dt = dt, T_max = T_max, use_cache = TRUE)
  n_total <- nrow(base$X)
  T_total <- base$T_max
  t_vec <- base$times

  # Define interior domain (0.05 T, 0.99 T) and boundary OOD domain [0, 0.05 T] U [0.99 T, T]
  idx_interior <- which(t_vec > 0.05 * T_total & t_vec < 0.99 * T_total)
  idx_ood_pool <- which(t_vec <= 0.05 * T_total | t_vec >= 0.99 * T_total)

  set.seed(as.integer(seed))

  # Sample train and non-overlapping in-sample test indices from interior
  n_needed_int <- min(length(idx_interior), N_train + N_test_in)
  sampled_int <- sample(idx_interior, size = n_needed_int, replace = FALSE)
  idx_train <- sort(sampled_int[seq_len(min(N_train, length(sampled_int)))])
  idx_test_in <- sort(sampled_int[(length(idx_train) + 1):length(sampled_int)])

  # Sample OOD test indices from boundary pool
  idx_ood_early <- which(t_vec <= 0.05 * T_total)
  idx_ood_late <- which(t_vec >= 0.99 * T_total)
  n_early <- min(length(idx_ood_early), floor(N_test_ood / 2))
  n_late <- min(length(idx_ood_late), N_test_ood - n_early)
  idx_test_ood <- sort(c(
    sample(idx_ood_early, size = n_early, replace = FALSE),
    sample(idx_ood_late, size = n_late, replace = FALSE)
  ))

  X_clean <- base$X
  Y_clean <- base$dX_fd
  m <- ncol(X_clean)

  if (noise_regime == "derivative_snr") {
    sd_deriv <- apply(Y_clean[idx_interior, , drop = FALSE], 2, sd)
    sd_deriv[sd_deriv < 1e-8] <- 1.0
    noise_sd_vec <- noise_eta * sd_deriv
    eps_train <- sapply(seq_len(m), function(j) rnorm(length(idx_train), 0, noise_sd_vec[j]))
    eps_in <- sapply(seq_len(m), function(j) rnorm(length(idx_test_in), 0, noise_sd_vec[j]))
    eps_ood <- sapply(seq_len(m), function(j) rnorm(length(idx_test_ood), 0, noise_sd_vec[j]))
    X_train <- X_clean[idx_train, , drop = FALSE]
    Y_train <- Y_clean[idx_train, , drop = FALSE] + eps_train
    X_test_in <- X_clean[idx_test_in, , drop = FALSE]
    Y_test_in_noisy <- Y_clean[idx_test_in, , drop = FALSE] + eps_in
    X_test_ood <- X_clean[idx_test_ood, , drop = FALSE]
    Y_test_ood_noisy <- Y_clean[idx_test_ood, , drop = FALSE] + eps_ood
  } else if (noise_regime == "iwsm_absolute") {
    noise_sd_vec <- rep(noise_eta, m)
    eps_train <- matrix(rnorm(length(idx_train) * m, 0, noise_eta), nrow = length(idx_train), ncol = m)
    eps_in <- matrix(rnorm(length(idx_test_in) * m, 0, noise_eta), nrow = length(idx_test_in), ncol = m)
    eps_ood <- matrix(rnorm(length(idx_test_ood) * m, 0, noise_eta), nrow = length(idx_test_ood), ncol = m)
    X_train <- X_clean[idx_train, , drop = FALSE]
    Y_train <- Y_clean[idx_train, , drop = FALSE] + eps_train
    X_test_in <- X_clean[idx_test_in, , drop = FALSE]
    Y_test_in_noisy <- Y_clean[idx_test_in, , drop = FALSE] + eps_in
    X_test_ood <- X_clean[idx_test_ood, , drop = FALSE]
    Y_test_ood_noisy <- Y_clean[idx_test_ood, , drop = FALSE] + eps_ood
  } else if (noise_regime == "state_noise") {
    # Subsample an equispaced trajectory segment of size 2500 from interior, add state noise, and differentiate
    idx_sub <- round(seq(min(idx_interior), max(idx_interior), length.out = max(N_train + N_test_in, 2000)))
    t_sub <- t_vec[idx_sub]
    X_sub_clean <- X_clean[idx_sub, , drop = FALSE]
    sd_state <- apply(X_sub_clean, 2, sd)
    sd_state[sd_state < 1e-8] <- 1.0
    noise_sd_vec <- noise_eta * sd_state
    nu_mat <- sapply(seq_len(m), function(j) rnorm(nrow(X_sub_clean), 0, noise_sd_vec[j]))
    X_sub_noisy <- X_sub_clean + nu_mat
    Y_sub_est <- estimate_derivatives_sgolay_spline(X_sub_noisy, t_sub)

    perm <- sample.int(nrow(X_sub_noisy))
    tr_loc <- sort(perm[seq_len(min(N_train, length(perm)))])
    in_loc <- sort(perm[(length(tr_loc) + 1):min(length(tr_loc) + N_test_in, length(perm))])

    X_train <- X_sub_noisy[tr_loc, , drop = FALSE]
    Y_train <- Y_sub_est[tr_loc, , drop = FALSE]
    idx_train <- idx_sub[tr_loc]

    X_test_in <- X_clean[idx_sub[in_loc], , drop = FALSE]
    Y_test_in_noisy <- Y_sub_est[in_loc, , drop = FALSE]
    idx_test_in <- idx_sub[in_loc]

    X_test_ood <- X_clean[idx_test_ood, , drop = FALSE]
    Y_test_ood_noisy <- Y_clean[idx_test_ood, , drop = FALSE]
  }

  # Rollout validation segment (continuous trajectory segment from test domain)
  rollout_start_idx <- min(idx_ood_late)
  rollout_indices <- seq(rollout_start_idx, min(n_total, rollout_start_idx + 5000), by = 10)

  list(
    system = base$system,
    seed = seed,
    N_train = nrow(X_train),
    noise_eta = noise_eta,
    noise_regime = noise_regime,
    noise_sd_vec = noise_sd_vec,
    X_train = X_train,
    Y_train = Y_train,
    Y_train_clean = Y_clean[idx_train, , drop = FALSE],
    X_test_in = X_test_in,
    Y_test_in_clean = Y_clean[idx_test_in, , drop = FALSE],
    Y_test_in_noisy = Y_test_in_noisy,
    X_test_ood = X_test_ood,
    Y_test_ood_clean = Y_clean[idx_test_ood, , drop = FALSE],
    Y_test_ood_noisy = Y_test_ood_noisy,
    rollout = list(
      times = t_vec[rollout_indices] - t_vec[rollout_indices[1]],
      X_true = X_clean[rollout_indices, , drop = FALSE],
      x0 = X_clean[rollout_indices[1], ]
    )
  )
}

