# ==============================================================================
# R/03_fbms_wrappers.R
# Tuned BGNLM (gmjmcmc.parallel) & Predefined-Basis BGLM (mjmcmc.parallel) Wrappers
# ==============================================================================

suppressPackageStartupMessages({
  library(FBMS)
  library(yaml)
})

if (!exists("build_poly_library", mode = "function")) {
  if (file.exists("R/02_basis_builders.R")) source("R/02_basis_builders.R")
  else if (file.exists("../../R/02_basis_builders.R")) source("../../R/02_basis_builders.R")
}
if (!exists("extract_canonical_atoms", mode = "function")) {
  if (file.exists("R/05_symbolic_canonicalize.R")) source("R/05_symbolic_canonicalize.R")
  else if (file.exists("../../R/05_symbolic_canonicalize.R")) source("../../R/05_symbolic_canonicalize.R")
}

#' Register radian-accurate physical operator transformations in .GlobalEnv for FBMS
#' and patch FBMS:::check.collinearity for low-dimensional ODE state spaces (m = 2, 3, 5)
register_fbms_transforms <- function() {
  assign("sin_rad", function(x) sin(x), envir = .GlobalEnv)
  assign("cos_rad", function(x) cos(x), envir = .GlobalEnv)
  assign("exp_neg", function(x) exp(-abs(x)), envir = .GlobalEnv)
  assign("inv_pos", function(x) 1.0 / (1.0 + abs(x)), envir = .GlobalEnv)
  assign("inv_sq",  function(x) 1.0 / (1.0 + x^2), envir = .GlobalEnv)
  assign("hill_2",  function(x) (x^2) / (1.0 + x^2), envir = .GlobalEnv)

  # Patch FBMS:::check.collinearity so n >= 3 * length(features) + 10 for low-D dynamical systems
  fix_col <- function(proposal, features, F.0.size, data, mock) {
    features[[length(features) + 1L]] <- proposal
    n <- min(nrow(data$x), max(F.0.size * 2L, length(features) * 3L + 15L))
    if (mock) {
      mock.data <- list(
        x = matrix(runif(n * (F.0.size + data$fixed), -100, 100), n, F.0.size + data$fixed),
        y = matrix(runif(n * ncol(as.matrix(data$y)), -100, 100), n, ncol(as.matrix(data$y))),
        fixed = data$fixed
      )
    } else {
      obs_idx <- seq_len(n)
      mock.data <- list(
        x = data$x[obs_idx, , drop = FALSE],
        y = as.matrix(data$y)[obs_idx, , drop = FALSE],
        fixed = data$fixed
      )
    }
    mock.data.precalc <- FBMS:::precalc.features(mock.data, features)
    linearmod <- lm.fit(mock.data.precalc$x, mock.data.precalc$y)
    if (sum(is.na(linearmod$coefficients)) == 0L) return(FALSE) else return(TRUE)
  }
  environment(fix_col) <- asNamespace("FBMS")
  tryCatch(assignInNamespace("check.collinearity", fix_col, ns = "FBMS"), error = function(e) NULL)

  invisible(list(
    transforms_core = c("sin_rad", "cos_rad", "p2", "p3", "inv_pos", "inv_sq"),
    transforms_extended = c("sin_rad", "cos_rad", "p2", "p3", "p05", "pm1", "exp_neg", "inv_pos", "inv_sq", "hill_2", "sigmoid"),
    transforms_iwsm = c("sin_deg", "cos_deg", "p0", "p2", "p3", "p05", "pm1", "pm2", "pm05")
  ))
}

#' Load Tuned BGNLM Hyperparameters from config/best_bgnlm_params.yaml (or return defaults)
get_tuned_bgnlm_config <- function(N_obs = 500L, p_cov = 3L) {
  cfg_paths <- c("config/best_bgnlm_params.yaml", "../../config/best_bgnlm_params.yaml")
  cfg_file <- cfg_paths[file.exists(cfg_paths)][1]

  if (!is.na(cfg_file) && file.exists(cfg_file)) {
    cfg <- yaml::read_yaml(cfg_file)
  } else {
    cfg <- list(
      beta_prior_type = "g-prior",
      g_rule = "max_n_p2",
      var = "unknown",
      r_rule = "N^(-1.0)",
      p_inter = 0.55,
      p_mod = 0.35,
      p_proj = 0.00,
      p_new = 0.10,
      pop_max = 20L,
      max_depth_D = 2L,
      P_generations = 15L,
      N_mcmc = 500L,
      runs = 6L,
      transforms_tier = "transforms_extended",
      energy_prune_tol = 0.025
    )
  }

  # Evaluate beta_prior
  if (cfg$beta_prior_type == "g-prior") {
    g_val <- max(as.numeric(N_obs), as.numeric(p_cov)^2)
    beta_prior <- list(type = "g-prior", g = g_val)
  } else if (cfg$beta_prior_type == "Jeffreys-BIC") {
    beta_prior <- list(type = "Jeffreys-BIC", var = cfg$var)
  } else {
    beta_prior <- list(type = cfg$beta_prior_type)
  }

  # Evaluate complexity penalty r
  N <- as.numeric(N_obs)
  r_val <- switch(
    as.character(cfg$r_rule),
    "N^(-0.5)" = N^(-0.5),
    "N^(-1.0)" = N^(-1.0),
    "N^(-1.5)" = N^(-1.5),
    "exp(-2)"  = exp(-2.0),
    "exp(-4)"  = exp(-4.0),
    1.0 / N
  )

  list(
    beta_prior = beta_prior,
    model_prior = list(r = r_val),
    p_inter = as.numeric(cfg$p_inter),
    p_mod = as.numeric(cfg$p_mod),
    p_proj = as.numeric(cfg$p_proj),
    p_new = as.numeric(cfg$p_new),
    pop_max = as.integer(cfg$pop_max),
    max_depth_D = as.integer(cfg$max_depth_D),
    P_generations = as.integer(cfg$P_generations),
    N_mcmc = as.integer(cfg$N_mcmc),
    runs = as.integer(cfg$runs),
    transforms_tier = as.character(cfg$transforms_tier),
    energy_prune_tol = if (!is.null(cfg$energy_prune_tol)) as.numeric(cfg$energy_prune_tol) else 0.025
  )
}

#' Evaluate an FBMS feature string safely on a state/design matrix
eval_fbms_feature_string <- function(feat_str, X_df) {
  register_fbms_transforms()
  s <- trimws(as.character(feat_str))
  # Normalize x.1 -> x1 if needed
  for (k in seq_len(ncol(X_df))) {
    s <- gsub(paste0("\\bx\\.", k, "\\b"), colnames(X_df)[k], s)
  }
  env <- list2env(as.list(X_df), parent = .GlobalEnv)
  val <- tryCatch(
    eval(parse(text = s), envir = env),
    error = function(e) {
      if (s %in% colnames(X_df)) return(X_df[[s]])
      evaluate_basis_atom(s, X_df)
    }
  )
  val <- as.numeric(val)
  val[!is.finite(val)] <- 0.0
  val
}

#' Fit a single FBMS method across all m equations of a dynamical system
run_fbms_method <- function(
  method_id = c(
    "BGNLM_Full_Tuned",
    "BGNLM_IWSM_Untuned",
    "BGNLM_Hybrid_Poly",
    "BGLM_Oracle_Basis",
    "BGLM_Poly_Basis",
    "BGLM_Overcomplete"
  ),
  task_data,
  cores = 2L,
  tuning_override = NULL,
  fast_mode = FALSE
) {
  method_id <- match.arg(method_id)
  t_dicts <- register_fbms_transforms()

  X_tr_raw <- as.matrix(task_data$X_train)
  Y_tr <- as.matrix(task_data$Y_train)
  X_in_raw <- as.matrix(task_data$X_test_in)
  X_ood_raw <- as.matrix(task_data$X_test_ood)
  m <- ncol(Y_tr)
  N_obs <- nrow(X_tr_raw)
  extra_oracle <- if (!is.null(task_data$system$oracle_extra_basis)) task_data$system$oracle_extra_basis else character(0)

  # Build design matrix X depending on method_id
  build_X_fn <- function(X_mat) {
    if (method_id %in% c("BGNLM_Full_Tuned", "BGNLM_IWSM_Untuned")) {
      return(X_mat)
    } else if (method_id == "BGNLM_Hybrid_Poly") {
      return(build_poly_library(X_mat, max_deg = 2L))
    } else if (method_id == "BGLM_Oracle_Basis") {
      return(build_oracle_library(X_mat, extra_atoms = extra_oracle, max_deg = 3L))
    } else if (method_id == "BGLM_Poly_Basis") {
      return(build_poly_library(X_mat, max_deg = 3L))
    } else if (method_id == "BGLM_Overcomplete") {
      return(build_overcomplete_library(X_mat, extra_atoms = extra_oracle))
    }
  }

  X_tr <- build_X_fn(X_tr_raw)
  X_in <- build_X_fn(X_in_raw)
  X_ood <- build_X_fn(X_ood_raw)
  orig_labels <- colnames(X_tr)

  # Safe R syntactic column names for FBMS formula interface, mapped back via orig_labels
  safe_cols <- paste0("v", seq_len(ncol(X_tr)))
  X_tr_df <- as.data.frame(X_tr); colnames(X_tr_df) <- safe_cols
  X_in_df <- as.data.frame(X_in); colnames(X_in_df) <- safe_cols
  X_ood_df <- as.data.frame(X_ood); colnames(X_ood_df) <- safe_cols

  cfg <- get_tuned_bgnlm_config(N_obs = N_obs, p_cov = ncol(X_tr))
  if (!is.null(tuning_override)) {
    for (nm in names(tuning_override)) cfg[[nm]] <- tuning_override[[nm]]
  }
  if (fast_mode) {
    cfg$P_generations <- min(cfg$P_generations, 6L)
    cfg$N_mcmc <- min(cfg$N_mcmc, 250L)
    cfg$runs <- min(cfg$runs, 2L)
  }

  discovered_mpm <- vector("list", m)
  discovered_map <- vector("list", m)
  pip_list <- vector("list", m)
  eq_models <- vector("list", m)

  pred_train <- matrix(0.0, nrow = nrow(X_tr), ncol = m)
  pred_in <- matrix(0.0, nrow = nrow(X_in), ncol = m)
  pred_ood <- matrix(0.0, nrow = nrow(X_ood), ncol = m)
  lower_in <- matrix(0.0, nrow = nrow(X_in), ncol = m)
  upper_in <- matrix(0.0, nrow = nrow(X_in), ncol = m)
  lower_ood <- matrix(0.0, nrow = nrow(X_ood), ncol = m)
  upper_ood <- matrix(0.0, nrow = nrow(X_ood), ncol = m)

  for (j in seq_len(m)) {
    df_j <- cbind(yc = as.numeric(Y_tr[, j]), X_tr_df)
    form_j <- as.formula(paste("yc ~", paste(safe_cols, collapse = " + ")))

    if (method_id == "BGNLM_IWSM_Untuned") {
      # Exact IWSM baseline settings (untuned)
      params <- FBMS::gen.params.gmjmcmc(ncol(X_tr_df))
      params$feat$pop.max <- 15L
      params$feat$keep.org <- TRUE
      params$feat$check.col <- TRUE
      P_use <- if (fast_mode) 5L else 12L
      N_use <- if (fast_mode) 200L else 400L
      runs_use <- if (fast_mode) 2L else 4L
      fbms_fit <- suppressWarnings(FBMS::fbms(
        formula = form_j,
        family = "gaussian",
        beta_prior = list(type = "Jeffreys-BIC", var = 1.0),
        model_prior = list(r = 1.0 / nrow(df_j)),
        data = df_j,
        method = "gmjmcmc.parallel",
        params = params,
        runs = runs_use,
        cores = max(1L, as.integer(cores)),
        N = N_use,
        P = P_use,
        transforms = t_dicts$transforms_iwsm,
        verbose = FALSE
      ))
      prune_tol <- 1e-4 # Untuned IWSM did not prune weak features
    } else if (method_id %in% c("BGNLM_Full_Tuned", "BGNLM_Hybrid_Poly")) {
      transforms_vec <- if (cfg$transforms_tier == "transforms_core") t_dicts$transforms_core else t_dicts$transforms_extended
      params <- FBMS::gen.params.gmjmcmc(ncol(X_tr_df))
      params$feat$pop.max <- cfg$pop_max
      params$feat$D <- cfg$max_depth_D
      params$feat$keep.org <- TRUE
      params$feat$check.col <- TRUE
      params$feat$eps <- 0.20
      probs <- FBMS::gen.probs.gmjmcmc(transforms_vec)
      probs$gen <- c(cfg$p_inter, cfg$p_mod, cfg$p_proj, cfg$p_new)

      use_custom_g <- identical(cfg$beta_prior$type, "g-prior")
      fbms_fit <- suppressWarnings(FBMS::fbms(
        formula = form_j,
        family = if (use_custom_g) "custom" else "gaussian",
        loglik.pi = if (use_custom_g) FBMS:::gaussian.loglik.g else NULL,
        beta_prior = if (use_custom_g) list(g = cfg$beta_prior$g) else cfg$beta_prior,
        model_prior = cfg$model_prior,
        data = df_j,
        method = "gmjmcmc.parallel",
        params = params,
        probs = probs,
        runs = cfg$runs,
        cores = max(1L, as.integer(cores)),
        N = cfg$N_mcmc,
        P = cfg$P_generations,
        transforms = transforms_vec,
        verbose = FALSE
      ))
      prune_tol <- cfg$energy_prune_tol
    } else {
      # Predefined-Basis BGLM via mjmcmc.parallel
      N_mcmc_bglm <- if (fast_mode) 400L else max(800L, cfg$N_mcmc)
      runs_bglm <- if (fast_mode) 2L else max(4L, cfg$runs)
      use_custom_g <- identical(cfg$beta_prior$type, "g-prior")
      fbms_fit <- suppressWarnings(FBMS::fbms(
        formula = form_j,
        family = if (use_custom_g) "custom" else "gaussian",
        loglik.pi = if (use_custom_g) FBMS:::gaussian.loglik.g else NULL,
        beta_prior = if (use_custom_g) list(g = cfg$beta_prior$g) else cfg$beta_prior,
        model_prior = cfg$model_prior,
        data = df_j,
        method = "mjmcmc.parallel",
        runs = runs_bglm,
        cores = max(1L, as.integer(cores)),
        N = N_mcmc_bglm,
        verbose = FALSE
      ))
      prune_tol <- cfg$energy_prune_tol
    }

    # Extract full PIP table with original column labels
    sum_all <- suppressWarnings(summary(fbms_fit, tol = 0.0, labels = orig_labels, verbose = FALSE))
    feat_strs_all <- as.character(sum_all$feats.strings)
    pips_all <- as.numeric(sum_all$marg.probs)
    names(pips_all) <- feat_strs_all
    pip_list[[j]] <- pips_all

    # Extract MAP model via get.best.model
    best_mod <- tryCatch(
      suppressWarnings(FBMS::get.best.model(fbms_fit, labels = orig_labels)),
      error = function(e) NULL
    )
    if (!is.null(best_mod) && length(best_mod$coefs) > 1L) {
      active_map <- names(best_mod$coefs)[-1][best_mod$coefs[-1] != 0]
      discovered_map[[j]] <- extract_canonical_atoms(active_map)
    } else {
      discovered_map[[j]] <- character(0)
    }

    # Extract MPM features (PIP > 0.5)
    mpm_idx <- which(pips_all > 0.5)
    mpm_raw_feats <- feat_strs_all[mpm_idx]

    # Evaluate selected MPM features on train/in/ood and refit posterior mode / OLS with energy pruning
    if (length(mpm_raw_feats) > 0L) {
      X_tr_named <- as.data.frame(X_tr); colnames(X_tr_named) <- orig_labels
      Phi_mpm_tr <- sapply(mpm_raw_feats, function(fs) eval_fbms_feature_string(fs, X_tr_named))
      if (!is.matrix(Phi_mpm_tr)) Phi_mpm_tr <- matrix(Phi_mpm_tr, nrow = nrow(X_tr_named))
      colnames(Phi_mpm_tr) <- mpm_raw_feats

      # Energy pruning of negligible MPM terms
      fit_init <- lm.fit(x = cbind(1.0, Phi_mpm_tr), y = Y_tr[, j])
      coefs_init <- fit_init$coefficients[-1]
      coefs_init[is.na(coefs_init)] <- 0.0
      y_norm <- max(sqrt(sum(Y_tr[, j]^2)), 1e-12)
      rel_energy <- sapply(seq_along(mpm_raw_feats), function(k) {
        abs(coefs_init[k]) * sqrt(sum(Phi_mpm_tr[, k]^2)) / y_norm
      })
      keep_k <- which(rel_energy >= prune_tol)
      mpm_kept_feats <- mpm_raw_feats[keep_k]
    } else {
      mpm_kept_feats <- character(0)
    }

    discovered_mpm[[j]] <- extract_canonical_atoms(mpm_kept_feats)

    # Fit final linear predictor on kept MPM features + compute BMA predictive intervals
    X_tr_named <- as.data.frame(X_tr); colnames(X_tr_named) <- orig_labels
    X_in_named <- as.data.frame(X_in); colnames(X_in_named) <- orig_labels
    X_ood_named <- as.data.frame(X_ood); colnames(X_ood_named) <- orig_labels

    if (length(mpm_kept_feats) > 0L) {
      Phi_tr <- as.matrix(sapply(mpm_kept_feats, function(fs) eval_fbms_feature_string(fs, X_tr_named)))
      Phi_in <- as.matrix(sapply(mpm_kept_feats, function(fs) eval_fbms_feature_string(fs, X_in_named)))
      Phi_ood <- as.matrix(sapply(mpm_kept_feats, function(fs) eval_fbms_feature_string(fs, X_ood_named)))
      if (nrow(X_tr_named) == 1L) Phi_tr <- matrix(Phi_tr, nrow = 1L)
      if (nrow(X_in_named) == 1L) Phi_in <- matrix(Phi_in, nrow = 1L)
      if (nrow(X_ood_named) == 1L) Phi_ood <- matrix(Phi_ood, nrow = 1L)

      X_des_tr <- cbind(1.0, Phi_tr)
      X_des_in <- cbind(1.0, Phi_in)
      X_des_ood <- cbind(1.0, Phi_ood)
      fit_final <- lm.fit(x = X_des_tr, y = Y_tr[, j])
      beta_hat <- fit_final$coefficients
      beta_hat[is.na(beta_hat)] <- 0.0
      df_res <- max(1L, nrow(X_des_tr) - qr(X_des_tr)$rank)
      sigma_hat <- sqrt(sum(fit_final$residuals^2) / df_res)

      mu_tr <- as.numeric(X_des_tr %*% beta_hat)
      mu_in <- as.numeric(X_des_in %*% beta_hat)
      mu_ood <- as.numeric(X_des_ood %*% beta_hat)

      # Epistemic + Aleatoric variance for 95% credible predictive intervals
      XtX_inv <- tryCatch(chol2inv(chol(crossprod(X_des_tr) + diag(1e-6, ncol(X_des_tr)))),
                          error = function(e) diag(1.0 / nrow(X_des_tr), ncol(X_des_tr)))
      se_in <- sigma_hat * sqrt(1.0 + pmax(0.0, rowSums((X_des_in %*% XtX_inv) * X_des_in)))
      se_ood <- sigma_hat * sqrt(1.0 + pmax(0.0, rowSums((X_des_ood %*% XtX_inv) * X_des_ood)))
    } else {
      beta_hat <- mean(Y_tr[, j])
      sigma_hat <- sd(Y_tr[, j])
      mu_tr <- rep(beta_hat, nrow(X_tr))
      mu_in <- rep(beta_hat, nrow(X_in))
      mu_ood <- rep(beta_hat, nrow(X_ood))
      se_in <- rep(sigma_hat, nrow(X_in))
      se_ood <- rep(sigma_hat, nrow(X_ood))
    }

    # Combine with FBMS BMA quantiles if available
    bma_in <- tryCatch(suppressWarnings(predict(fbms_fit, x = X_in_df)), error = function(e) NULL)
    bma_ood <- tryCatch(suppressWarnings(predict(fbms_fit, x = X_ood_df)), error = function(e) NULL)

    if (!is.null(bma_in) && !is.null(bma_in$aggr$mean) && all(is.finite(bma_in$aggr$mean))) {
      # Blend BMA mean and MPM mean, and widen BMA quantiles by residual dispersion sigma_hat
      q_low_bma <- bma_in$aggr$quantiles[1, ]
      q_high_bma <- bma_in$aggr$quantiles[3, ]
      bma_halfwidth <- 0.5 * (q_high_bma - q_low_bma)
      tot_halfwidth_in <- sqrt(bma_halfwidth^2 + (1.96 * se_in)^2)
      lower_in[, j] <- mu_in - tot_halfwidth_in
      upper_in[, j] <- mu_in + tot_halfwidth_in
    } else {
      lower_in[, j] <- mu_in - 1.96 * se_in
      upper_in[, j] <- mu_in + 1.96 * se_in
    }

    if (!is.null(bma_ood) && !is.null(bma_ood$aggr$mean) && all(is.finite(bma_ood$aggr$mean))) {
      q_low_bma <- bma_ood$aggr$quantiles[1, ]
      q_high_bma <- bma_ood$aggr$quantiles[3, ]
      bma_halfwidth <- 0.5 * (q_high_bma - q_low_bma)
      tot_halfwidth_ood <- sqrt(bma_halfwidth^2 + (1.96 * se_ood)^2)
      lower_ood[, j] <- mu_ood - tot_halfwidth_ood
      upper_ood[, j] <- mu_ood + tot_halfwidth_ood
    } else {
      lower_ood[, j] <- mu_ood - 1.96 * se_ood
      upper_ood[, j] <- mu_ood + 1.96 * se_ood
    }

    pred_train[, j] <- mu_tr
    pred_in[, j] <- mu_in
    pred_ood[, j] <- mu_ood

    eq_models[[j]] <- list(
      feats = mpm_kept_feats,
      beta = beta_hat,
      sigma = sigma_hat
    )
  }

  # Construct fast RHS predictor for ODE trajectory rollout
  predict_rhs_fn <- function(X_state_mat) {
    X_full <- build_X_fn(X_state_mat)
    X_df_new <- as.data.frame(X_full)
    colnames(X_df_new) <- orig_labels
    out_d <- matrix(0.0, nrow = nrow(X_df_new), ncol = m)
    for (j in seq_len(m)) {
      mod_j <- eq_models[[j]]
      if (length(mod_j$feats) == 0L) {
        out_d[, j] <- mod_j$beta[1]
      } else {
        Phi_new <- sapply(mod_j$feats, function(fs) eval_fbms_feature_string(fs, X_df_new))
        if (nrow(X_df_new) == 1L) Phi_new <- matrix(Phi_new, nrow = 1L)
        out_d[, j] <- as.numeric(cbind(1.0, Phi_new) %*% mod_j$beta)
      }
    }
    out_d
  }

  list(
    method_id = method_id,
    discovered = discovered_mpm,
    discovered_map = discovered_map,
    pip_list = pip_list,
    pred_train = pred_train,
    pred_in = pred_in,
    pred_ood = pred_ood,
    lower_in = lower_in,
    upper_in = upper_in,
    lower_ood = lower_ood,
    upper_ood = upper_ood,
    predict_rhs_fn = predict_rhs_fn
  )
}

