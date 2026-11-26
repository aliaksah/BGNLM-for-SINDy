library(testthat)
source("../../R/01_ode_simulators.R")
source("../../R/02_basis_builders.R")
source("../../R/05_symbolic_canonicalize.R")
source("../../R/06_metrics_consistency.R")

test_that("Algebraically equivalent strings across FBMS, PySR, PyOperon, and PySINDy canonicalize identically", {
  # Commutative products
  expect_equal(extract_canonical_atoms("x2 * x1"), "x1*x2")
  expect_equal(extract_canonical_atoms("(x.2 * x.1)"), "x1*x2")
  expect_equal(extract_canonical_atoms("x1 x2"), "x1*x2")

  # Powers and FBMS transforms
  expect_equal(extract_canonical_atoms("p2(x1) * x2"), "x1^2*x2")
  expect_equal(extract_canonical_atoms("x2 * (x1^2)"), "x1^2*x2")
  expect_equal(extract_canonical_atoms("square(x_1) * x_2"), "x1^2*x2")
  expect_equal(extract_canonical_atoms("x1**2 * x2"), "x1^2*x2")
  expect_equal(extract_canonical_atoms("p3(x.1)"), "x1^3")
  expect_equal(extract_canonical_atoms("cube(x1)"), "x1^3")

  # Trigonometric and fractional transforms (sin_deg is degree-scaled and MUST NOT match radian sin)
  expect_equal(extract_canonical_atoms("sin_rad(x1)"), "sin(x1)")
  expect_false(identical(extract_canonical_atoms("sin_deg(x.1)"), "sin(x1)"))
  expect_equal(extract_canonical_atoms("p05(x.1)"), "sqrt(x1)")
  expect_equal(extract_canonical_atoms("sqroot(x1)"), "sqrt(x1)")
  expect_equal(extract_canonical_atoms("x1^(-1)"), "1/(x1)")
  expect_equal(extract_canonical_atoms("x1^(0.5)"), "sqrt(x1)")

  # Rational / Hill transforms and rational complement equivalence modulo intercept
  expect_equal(extract_canonical_atoms("inv_sq(x.3)"), "1/(1+x3^2)")
  expect_equal(extract_canonical_atoms("x.3 * inv_pos(x.3)"), "x3/(1+x3)")
  expect_equal(extract_canonical_atoms("x.1 * hill_2(x.2)"), "x1*x2^2/(1+x2^2)")
  expect_equal(normalize_modulo_intercept("x2^2/(1+x2^2)"), "1/(1+x2^2)")
  expect_equal(normalize_modulo_intercept("x1/(1+x1)"), "1/(1+x1)")

  # Full expanded equation from PySR / PyOperon
  expect_equal(
    extract_canonical_atoms("x1 * (28.0 - x3) - 1.0 * x2"),
    c("x1", "x1*x3", "x2")
  )
})

test_that("Selection consistency metrics produce TPR=1, FDR=0, Exact_Recovery=1 for equivalent representations", {
  gt_lorenz <- list(
    c("x1", "x2"),
    c("x1", "x2", "x1*x3"),
    c("x3", "x1*x2")
  )
  discovered_fbms <- list(
    c("x.2", "x.1"),
    c("(x.3 * x.1)", "x.1", "x.2"),
    c("(x.2*x.1)", "x.3")
  )
  m <- compute_structure_metrics(discovered_fbms, gt_lorenz)
  expect_equal(m$TPR, 1.0)
  expect_equal(m$FDR, 0.0)
  expect_equal(m$F1, 1.0)
  expect_equal(m$Exact_Recovery, 1L)
  expect_equal(m$SHD, 0L)
})

test_that("Jaccard and Nogueira stability estimators work properly", {
  sets_perfect <- replicate(10, c("eq1::x1", "eq1::x2", "eq2::x1*x3"), simplify = FALSE)
  expect_equal(compute_jaccard_stability(sets_perfect), 1.0)
  nog <- compute_nogueira_stability(sets_perfect, d_total = 100)
  expect_equal(nog$phi, 1.0)
})

test_that("STLSQ returns lambda and PIP calibration metrics are well-formed", {
  source("../../R/04_sindy_baselines.R")
  set.seed(1)
  X <- matrix(rnorm(100 * 2), ncol = 2); colnames(X) <- c("x1", "x2")
  Phi <- build_poly_library(X, max_deg = 2L)
  y <- 2.0 * X[, 1] - 1.5 * X[, 1] * X[, 2] + rnorm(100, 0, 0.05)
  cv_fit <- cv_stlsq_fit(Phi, y)
  expect_true(!is.null(cv_fit$lambda) && is.finite(cv_fit$lambda))

  pips <- list(c("x1" = 0.99, "x2" = 0.02, "x1*x2" = 0.98))
  gt <- list(c("x1", "x1*x2"))
  cal <- compute_pip_calibration_metrics(pips, gt)
  expect_true(is.finite(cal$PIP_SSPE) && cal$PIP_SSPE < 0.01)
  expect_equal(cal$PIP_AUROC, 1.0)
})

