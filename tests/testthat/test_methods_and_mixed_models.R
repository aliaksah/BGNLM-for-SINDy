library(testthat)
source("../../R/01_ode_simulators.R")
source("../../R/02_basis_builders.R")
source("../../R/05_symbolic_canonicalize.R")
source("../../R/06_metrics_consistency.R")
source("../../R/03_fbms_wrappers.R")
source("../../R/04_sindy_baselines.R")
source("../../R/04b_sr_baselines.R")
source("../../R/08_mixed_effects_eval.R")

test_that("FBMS, SINDy, and Symbolic Regression wrappers return standardized output objects", {
  td <- generate_synthetic_task_data(
    system_id = "sys04_lotka_volterra",
    N_train = 150,
    noise_eta = 0.01,
    seed = 42,
    T_max = 10.0
  )

  res_bgnlm <- run_fbms_method("BGNLM_Full_Tuned", td, cores = 1L, fast_mode = TRUE)
  expect_equal(length(res_bgnlm$discovered), 2L)
  expect_true(is.function(res_bgnlm$predict_rhs_fn))
  st_bgnlm <- compute_structure_metrics(res_bgnlm$discovered, td$system$ground_truth)
  expect_gte(st_bgnlm$TPR, 0.5)

  res_sindy <- run_sindy_baseline("SINDy_STLSQ_Oracle", td)
  expect_equal(length(res_sindy$discovered), 2L)
  st_sindy <- compute_structure_metrics(res_sindy$discovered, td$system$ground_truth)
  expect_equal(st_sindy$TPR, 1.0)

  res_sr <- run_sr_baseline("SR_PySR", td)
  expect_equal(length(res_sr$discovered), 2L)
  expect_equal(dim(res_sr$pred_ood), dim(td$Y_test_ood_clean))
})

