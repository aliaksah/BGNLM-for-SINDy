# ==============================================================================
# R/10_real_data_experiments.R
# Real-World Empirical Datasets Preprocessing, Derivative Smoothing, & Temporal CV
# ==============================================================================

if (!exists("get_systems_registry", mode = "function")) {
  if (file.exists("R/01_ode_simulators.R")) {
    source("R/01_ode_simulators.R")
  } else if (file.exists("../../R/01_ode_simulators.R")) {
    source("../../R/01_ode_simulators.R")
  }
}

#' Registry of the 5 Canonical Real-World Empirical Datasets
get_real_datasets_registry <- function() {
  list(
    real01_lynx_hare = list(
      dataset_id = 1L,
      id = "real01_lynx_hare",
      name = "Hudson Bay Lynx-Hare Pelts",
      domain = "Ecology",
      file = "data/real_world/lynx_hare_hudson_bay.csv",
      dim = 2L,
      state_names = c("x1", "x2"),
      input_names = character(0),
      feature_names = c("x1", "x2"),
      ground_truth = list(
        c("x1", "x1*x2"),
        c("x2", "x1*x2")
      ),
      oracle_extra_basis = character(0)
    ),

    real02_paramecium_didinium = list(
      dataset_id = 2L,
      id = "real02_paramecium_didinium",
      name = "Gause-Veilleux Paramecium-Didinium",
      domain = "Experimental Microcosm",
      file = "data/real_world/gause_paramecium_didinium.csv",
      dim = 2L,
      state_names = c("x1", "x2"),
      input_names = character(0),
      feature_names = c("x1", "x2"),
      ground_truth = list(
        c("x1", "x1^2", "x1*x2/(1+x1)"),
        c("x2", "x1*x2/(1+x1)")
      ),
      oracle_extra_basis = c("1/(1+x1)", "x1/(1+x1)", "x1*x2/(1+x1)")
    ),

    real03_schmidt_pendulum = list(
      dataset_id = 3L,
      id = "real03_schmidt_pendulum",
      name = "Schmidt-Lipson Physical Pendulum",
      domain = "Lab Mechanics",
      file = "data/real_world/schmidt_real_pendulum.csv",
      dim = 2L,
      state_names = c("x1", "x2"),
      input_names = character(0),
      feature_names = c("x1", "x2"),
      ground_truth = list(
        c("x2"),
        c("x2", "sin(x1)")
      ),
      oracle_extra_basis = c("sin(x1)", "cos(x1)")
    ),

    real04_cascaded_tanks = list(
      dataset_id = 4L,
      id = "real04_cascaded_tanks",
      name = "Cascaded Water Tanks Benchmark",
      domain = "Nonlinear System Identification",
      file = "data/real_world/cascaded_tanks_benchmark.csv",
      dim = 2L,
      state_names = c("x1", "x2"),
      input_names = c("x3"),
      feature_names = c("x1", "x2", "x3"),
      ground_truth = list(
        c("sqrt(x1)", "x3"),
        c("sqrt(x1)", "sqrt(x2)")
      ),
      oracle_extra_basis = c("sqrt(x1)", "sqrt(x2)")
    ),

    real05_yeast_glycolysis = list(
      dataset_id = 5L,
      id = "real05_yeast_glycolysis",
      name = "Yeast Glycolytic NADH Oscillations",
      domain = "Biochemistry",
      file = "data/real_world/yeast_glycolysis_nadh.csv",
      dim = 2L,
      state_names = c("x1", "x2"),
      input_names = character(0),
      feature_names = c("x1", "x2"),
      ground_truth = list(
        c("x1", "x1*x2^2/(1+x2^2)"),
        c("x1", "x2", "x1*x2^2/(1+x2^2)")
      ),
      oracle_extra_basis = c("1/(1+x2^2)", "x2^2/(1+x2^2)", "x1*x2^2/(1+x2^2)")
    )
  )
}

#' Generate a Real-World Dataset Task Split (Rolling Temporal Split + 85% Subsampling)
generate_real_task_data <- function(
  dataset_id,
  diff_method = c("sgolay_spline", "tvregdiff"),
  fold_id = 1L,
  seed = 50101L
) {
  diff_method <- match.arg(diff_method)
  reg <- get_real_datasets_registry()
  if (!dataset_id %in% names(reg)) {
    stop("Unknown real dataset_id: ", dataset_id)
  }
  ds <- reg[[dataset_id]]
  if (!file.exists(ds$file)) {
    source("scripts/download_real_data.R")
  }
  df <- read.csv(ds$file, stringsAsFactors = FALSE)
  times <- as.numeric(df$time)
  X_states <- as.matrix(df[, ds$state_names, drop = FALSE])
  X_features <- as.matrix(df[, ds$feature_names, drop = FALSE])

  # Estimate derivatives on full series using specified differentiation scheme
  if (diff_method == "sgolay_spline") {
    Y_est <- estimate_derivatives_sgolay_spline(X_states, times)
  } else {
    Y_est <- estimate_derivatives_tvregdiff(X_states, times)
  }

  n_total <- nrow(X_features)
  # Rolling expanding-window cutoff around 70% train (varies between 65% and 75% across folds)
  set.seed(as.integer(seed))
  frac_train <- 0.65 + 0.10 * ((fold_id - 1L) %% 10L) / 9.0
  n_train_window <- max(20L, floor(frac_train * n_total))
  idx_train_pool <- seq_len(n_train_window)
  idx_ood <- seq(from = n_train_window + 1L, to = n_total)

  # 85% subsample of training window for structural stability; remaining 15% serves as in-sample validation
  n_sub <- max(15L, floor(0.85 * length(idx_train_pool)))
  idx_train <- sort(sample(idx_train_pool, size = n_sub, replace = FALSE))
  idx_in <- setdiff(idx_train_pool, idx_train)
  if (length(idx_in) < 3L) {
    idx_in <- idx_train
  }

  # Rollout horizon on held-out OOD test window
  rollout_idx <- idx_ood[seq_len(min(length(idx_ood), 100L))]
  rollout_times <- times[rollout_idx] - times[rollout_idx[1]]

  list(
    system = list(
      sys_id = ds$dataset_id,
      id = ds$id,
      name = ds$name,
      group = ds$domain,
      dim = ds$dim,
      state_names = ds$state_names,
      feature_names = ds$feature_names,
      input_names = ds$input_names,
      ground_truth = ds$ground_truth,
      oracle_extra_basis = ds$oracle_extra_basis
    ),
    dataset_id = ds$id,
    diff_method = diff_method,
    fold_id = fold_id,
    seed = seed,
    N_train = length(idx_train),
    X_train = X_features[idx_train, , drop = FALSE],
    Y_train = Y_est[idx_train, , drop = FALSE],
    Y_train_clean = Y_est[idx_train, , drop = FALSE],
    X_test_in = X_features[idx_in, , drop = FALSE],
    Y_test_in_clean = Y_est[idx_in, , drop = FALSE],
    Y_test_in_noisy = Y_est[idx_in, , drop = FALSE],
    X_test_ood = X_features[idx_ood, , drop = FALSE],
    Y_test_ood_clean = Y_est[idx_ood, , drop = FALSE],
    Y_test_ood_noisy = Y_est[idx_ood, , drop = FALSE],
    rollout = list(
      times = rollout_times,
      X_true = X_states[rollout_idx, , drop = FALSE],
      X_inputs = if (length(ds$input_names) > 0) X_features[rollout_idx, ds$input_names, drop = FALSE] else NULL,
      x0 = X_states[rollout_idx[1], ]
    )
  )
}
