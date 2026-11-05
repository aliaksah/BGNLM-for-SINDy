library(testthat)
source("../../R/01_ode_simulators.R")
source("../../R/02_basis_builders.R")
source("../../R/10_real_data_experiments.R")

test_that("All 12 synthetic ODE systems integrate cleanly and match analytical derivatives with R2 > 0.9999", {
  reg <- get_systems_registry()
  expect_equal(length(reg), 12L)

  for (sys_id in names(reg)) {
    traj <- simulate_base_trajectory(sys_id, T_max = 5.0, use_cache = FALSE)
    expect_true(all(is.finite(traj$X)), info = paste("Non-finite states in", sys_id))
    expect_true(all(is.finite(traj$dX_analytical)), info = paste("Non-finite analytical deriv in", sys_id))

    # Verify finite difference vs analytical derivative has R2 > 0.9999 at dt = 1e-4
    for (j in seq_len(traj$system$dim)) {
      r2_fd <- cor(traj$dX_fd[, j], traj$dX_analytical[, j])^2
      expect_gt(r2_fd, 0.9999, label = paste(sys_id, "eq", j, "FD vs analytical R2"))

      # Verify ground-truth basis atoms in Oracle library yield R2 > 0.9999 against analytical derivative
      Phi_oracle <- build_oracle_library(traj$X[1:1000, , drop = FALSE], extra_atoms = traj$system$oracle_extra_basis)
      gt_atoms <- traj$system$ground_truth[[j]]
      expect_true(all(gt_atoms %in% colnames(Phi_oracle)),
                  info = paste(sys_id, "eq", j, "missing atoms in oracle:", paste(setdiff(gt_atoms, colnames(Phi_oracle)), collapse = ",")))
      fit_lm <- lm(traj$dX_analytical[1:1000, j] ~ Phi_oracle[, gt_atoms, drop = FALSE])
      r2_gt <- suppressWarnings(summary(fit_lm)$r.squared)
      expect_gt(r2_gt, 0.9999, label = paste(sys_id, "eq", j, "ground-truth basis R2"))
    }
  }
})

test_that("Overcomplete library has > 60 candidate columns and finite values", {
  traj <- simulate_base_trajectory("sys02_lorenz63", T_max = 2.0, use_cache = FALSE)
  Phi_over <- build_overcomplete_library(traj$X[1:200, , drop = FALSE])
  expect_gt(ncol(Phi_over), 60L)
  expect_true(all(is.finite(Phi_over)))
})
