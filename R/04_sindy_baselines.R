# ==============================================================================
# R/04_sindy_baselines.R
# Classical STLSQ SINDy, Bayesian/Ensemble UQ-SINDy, and Rational SINDy-PI
# ==============================================================================

if (!exists("build_poly_library", mode = "function")) {
  if (file.exists("R/02_basis_builders.R")) source("R/02_basis_builders.R")
  else if (file.exists("../../R/02_basis_builders.R")) source("../../R/02_basis_builders.R")
}
if (!exists("extract_canonical_atoms", mode = "function")) {
  if (file.exists("R/05_symbolic_canonicalize.R")) source("R/05_symbolic_canonicalize.R")
  else if (file.exists("../../R/05_symbolic_canonicalize.R")) source("../../R/05_symbolic_canonicalize.R")
}

#' Fit Sequential Thresholded Least Squares (STLSQ) on standardized design matrix Phi
stlsq_fit_single <- function(Phi, y, lambda = 0.05, alpha = 1e-5, max_iter = 15L) {
  n <- nrow(Phi)
  p <- ncol(Phi)
  y_mean <- mean(y)
  y_c <- y - y_mean

  col_means <- colMeans(Phi)
  col_sds <- apply(Phi, 2, sd)
  valid_cols <- which(is.finite(col_sds) & col_sds > 1e-10)
  if (length(valid_cols) == 0L) {
    return(list(intercept = y_mean, coefs = setNames(rep(0.0, p), colnames(Phi)), active = integer(0)))
  }

  Phi_std <- scale(Phi[, valid_cols, drop = FALSE], center = col_means[valid_cols], scale = col_sds[valid_cols])
  y_sd <- sd(y_c)
  if (!is.finite(y_sd) || y_sd < 1e-10) y_sd <- 1.0
  y_std <- y_c / y_sd

  active_loc <- seq_along(valid_cols)
  w_std <- rep(0.0, length(valid_cols))

  for (it in seq_len(max_iter)) {
    if (length(active_loc) == 0L) break
    X_a <- Phi_std[, active_loc, drop = FALSE]
    XtX <- crossprod(X_a) + diag(alpha * n, ncol(X_a))
    Xty <- crossprod(X_a, y_std)
    w_a <- tryCatch(as.numeric(solve(XtX, Xty)), error = function(e) rep(0.0, ncol(X_a)))
    w_std <- rep(0.0, length(valid_cols))
    w_std[active_loc] <- w_a
    new_active <- which(abs(w_std) >= lambda)
    if (identical(new_active, active_loc)) break
    active_loc <- new_active
  }

  # Unbiased OLS refit on active support
  coefs_full <- rep(0.0, p)
  names(coefs_full) <- colnames(Phi)
  active_orig <- valid_cols[active_loc]
  if (length(active_orig) > 0L) {
    X_act <- cbind(1.0, Phi[, active_orig, drop = FALSE])
    fit_ols <- lm.fit(x = X_act, y = y)
    b_ols <- fit_ols$coefficients
    b_ols[is.na(b_ols)] <- 0.0
    intercept <- b_ols[1]
    coefs_full[active_orig] <- b_ols[-1]
  } else {
    intercept <- y_mean
  }

  list(intercept = intercept, coefs = coefs_full, active = active_orig)
}

#' Cross-Validated + BIC-Regularized STLSQ Selection
cv_stlsq_fit <- function(Phi, y, lambda_grid = exp(seq(log(0.005), log(0.60), length.out = 20L)), k_folds = 5L) {
  n <- nrow(Phi)
  fold_ids <- rep(seq_len(k_folds), length.out = n)
  cv_scores <- numeric(length(lambda_grid))

  for (idx_l in seq_along(lambda_grid)) {
    lam <- lambda_grid[idx_l]
    mse_folds <- numeric(k_folds)
    k_terms <- 0.0
    for (f in seq_len(k_folds)) {
      tr <- which(fold_ids != f)
      va <- which(fold_ids == f)
      fit_f <- stlsq_fit_single(Phi[tr, , drop = FALSE], y[tr], lambda = lam)
      pred_va <- fit_f$intercept + as.numeric(Phi[va, , drop = FALSE] %*% fit_f$coefs)
      mse_folds[f] <- mean((y[va] - pred_va)^2)
      k_terms <- k_terms + length(fit_f$active) / k_folds
    }
    # Combine CV MSE with BIC sparsity penalty to prevent overfitting under noise
    cv_mse <- mean(mse_folds)
    cv_scores[idx_l] <- n * log(max(cv_mse, 1e-12)) + log(n) * k_terms
  }

  best_lam <- lambda_grid[which.min(cv_scores)]
  stlsq_fit_single(Phi, y, lambda = best_lam)
}

#' Build Rational / Implicit Dictionary for SINDy-PI (Kaheman et al. 2020)
build_sindy_pi_rational_library <- function(X) {
  X <- as.matrix(X)
  col_names <- colnames(X)
  if (is.null(col_names)) {
    col_names <- paste0("x", seq_len(ncol(X)))
    colnames(X) <- col_names
  }
  m <- ncol(X)
  Phi_poly <- build_poly_library(X, max_deg = 2L)

  rat_list <- list()
  rat_names <- character(0)
  add_rat <- function(vec, nm) {
    if (!nm %in% c(colnames(Phi_poly), rat_names)) {
      vec[!is.finite(vec)] <- 0.0
      rat_list[[length(rat_list) + 1L]] <<- as.numeric(vec)
      rat_names <<- c(rat_names, nm)
    }
  }

  for (i in seq_len(m)) {
    xi <- X[, i]; ni <- col_names[i]
    add_rat(1.0 / (1.0 + abs(xi)), paste0("1/(1+", ni, ")"))
    add_rat(xi / (1.0 + abs(xi)), paste0(ni, "/(1+", ni, ")"))
    add_rat(1.0 / (1.0 + xi^2), paste0("1/(1+", ni, "^2)"))
    add_rat((xi^2) / (1.0 + xi^2), paste0(ni, "^2/(1+", ni, "^2)"))
    for (j in seq_len(m)) {
      xj <- X[, j]; nj <- col_names[j]
      if (i != j) {
        add_rat((xi * xj) / (1.0 + abs(xi)), paste0(ni, "*", nj, "/(1+", ni, ")"))
        add_rat((xi * (xj^2)) / (1.0 + xj^2), paste0(ni, "*", nj, "^2/(1+", nj, "^2)"))
      }
    }
  }

  Phi_rat <- do.call(cbind, rat_list)
  colnames(Phi_rat) <- rat_names
  cbind(Phi_poly, Phi_rat)
}

#' Unified Runner for the 5 SINDy Baselines
run_sindy_baseline <- function(
  method_id = c(
    "SINDy_STLSQ_Oracle",
    "SINDy_STLSQ_Poly",
    "UQ_SINDy_Oracle",
    "UQ_SINDy_Poly",
    "SINDy_PI_Rational"
  ),
  task_data,
  n_boot = 80L
) {
  method_id <- match.arg(method_id)
  X_tr_raw <- as.matrix(task_data$X_train)
  Y_tr <- as.matrix(task_data$Y_train)
  X_in_raw <- as.matrix(task_data$X_test_in)
  X_ood_raw <- as.matrix(task_data$X_test_ood)
  m <- ncol(Y_tr)
  extra_oracle <- if (!is.null(task_data$system$oracle_extra_basis)) task_data$system$oracle_extra_basis else character(0)

  build_Phi_fn <- function(X_mat) {
    if (method_id %in% c("SINDy_STLSQ_Oracle", "UQ_SINDy_Oracle")) {
      build_oracle_library(X_mat, extra_atoms = extra_oracle, max_deg = 3L)
    } else if (method_id %in% c("SINDy_STLSQ_Poly", "UQ_SINDy_Poly")) {
      build_poly_library(X_mat, max_deg = 3L)
    } else if (method_id == "SINDy_PI_Rational") {
      build_sindy_pi_rational_library(X_mat)
    }
  }

  Phi_tr <- build_Phi_fn(X_tr_raw)
  Phi_in <- build_Phi_fn(X_in_raw)
  Phi_ood <- build_Phi_fn(X_ood_raw)
  cand_names <- colnames(Phi_tr)
  p <- ncol(Phi_tr)

  discovered <- vector("list", m)
  pip_list <- vector("list", m)
  eq_models <- vector("list", m)

  pred_train <- matrix(0.0, nrow = nrow(Phi_tr), ncol = m)
  pred_in <- matrix(0.0, nrow = nrow(Phi_in), ncol = m)
  pred_ood <- matrix(0.0, nrow = nrow(Phi_ood), ncol = m)
  lower_in <- matrix(0.0, nrow = nrow(Phi_in), ncol = m)
  upper_in <- matrix(0.0, nrow = nrow(Phi_in), ncol = m)
  lower_ood <- matrix(0.0, nrow = nrow(Phi_ood), ncol = m)
  upper_ood <- matrix(0.0, nrow = nrow(Phi_ood), ncol = m)

  for (j in seq_len(m)) {
    yj <- Y_tr[, j]

    if (method_id %in% c("SINDy_STLSQ_Oracle", "SINDy_STLSQ_Poly", "SINDy_PI_Rational")) {
      fit_j <- cv_stlsq_fit(Phi_tr, yj)
      active_names <- cand_names[fit_j$active]
      discovered[[j]] <- extract_canonical_atoms(active_names)
      pip_vec <- as.numeric(cand_names %in% active_names)
      names(pip_vec) <- cand_names
      pip_list[[j]] <- pip_vec

      pred_train[, j] <- fit_j$intercept + as.numeric(Phi_tr %*% fit_j$coefs)
      pred_in[, j] <- fit_j$intercept + as.numeric(Phi_in %*% fit_j$coefs)
      pred_ood[, j] <- fit_j$intercept + as.numeric(Phi_ood %*% fit_j$coefs)

      # Point-estimate methods have no calibrated Bayesian credible intervals (NA or residual plug-in)
      lower_in[, j] <- NA_real_; upper_in[, j] <- NA_real_
      lower_ood[, j] <- NA_real_; upper_ood[, j] <- NA_real_
      eq_models[[j]] <- list(intercept = fit_j$intercept, coefs = fit_j$coefs)

    } else {
      # UQ_SINDy (Ensemble / Spike-and-Slab Bagging SINDy: Fasel et al. 2022 / Hirsh et al. 2022)
      fit_base <- cv_stlsq_fit(Phi_tr, yj)
      n_obs <- nrow(Phi_tr)
      boot_coefs <- matrix(0.0, nrow = n_boot, ncol = p)
      boot_inter <- numeric(n_boot)
      boot_active <- matrix(0.0, nrow = n_boot, ncol = p)

      # Use CV-selected threshold scale across bootstrap replicates
      lam_base <- max(0.03, median(abs(fit_base$coefs[fit_base$active]), na.rm = TRUE) * 0.25)
      if (!is.finite(lam_base)) lam_base <- 0.05

      for (b in seq_len(n_boot)) {
        idx_b <- sample.int(n_obs, size = floor(0.85 * n_obs), replace = TRUE)
        fb <- stlsq_fit_single(Phi_tr[idx_b, , drop = FALSE], yj[idx_b], lambda = lam_base)
        boot_coefs[b, ] <- fb$coefs
        boot_inter[b] <- fb$intercept
        boot_active[b, fb$active] <- 1.0
      }

      pips <- colMeans(boot_active)
      names(pips) <- cand_names
      pip_list[[j]] <- pips

      # Median Probability Model (PIP > 0.5)
      mpm_idx <- which(pips > 0.5)
      discovered[[j]] <- extract_canonical_atoms(cand_names[mpm_idx])

      if (length(mpm_idx) > 0L) {
        X_act <- cbind(1.0, Phi_tr[, mpm_idx, drop = FALSE])
        ols <- lm.fit(x = X_act, y = yj)
        b_ols <- ols$coefficients; b_ols[is.na(b_ols)] <- 0.0
        inter_mpm <- b_ols[1]
        coefs_mpm <- rep(0.0, p); names(coefs_mpm) <- cand_names
        coefs_mpm[mpm_idx] <- b_ols[-1]
        sigma_res <- sqrt(mean(ols$residuals^2))
      } else {
        inter_mpm <- mean(yj)
        coefs_mpm <- setNames(rep(0.0, p), cand_names)
        sigma_res <- sd(yj)
      }

      pred_train[, j] <- inter_mpm + as.numeric(Phi_tr %*% coefs_mpm)
      pred_in[, j] <- inter_mpm + as.numeric(Phi_in %*% coefs_mpm)
      pred_ood[, j] <- inter_mpm + as.numeric(Phi_ood %*% coefs_mpm)

      # Ensemble predictive quantiles + residual standard deviation
      ens_in <- sweep(Phi_in %*% t(boot_coefs), 2, boot_inter, "+")
      ens_ood <- sweep(Phi_ood %*% t(boot_coefs), 2, boot_inter, "+")
      sd_ens_in <- apply(ens_in, 1, sd)
      sd_ens_ood <- apply(ens_ood, 1, sd)
      tot_sd_in <- sqrt(sd_ens_in^2 + sigma_res^2)
      tot_sd_ood <- sqrt(sd_ens_ood^2 + sigma_res^2)

      lower_in[, j] <- pred_in[, j] - 1.96 * tot_sd_in
      upper_in[, j] <- pred_in[, j] + 1.96 * tot_sd_in
      lower_ood[, j] <- pred_ood[, j] - 1.96 * tot_sd_ood
      upper_ood[, j] <- pred_ood[, j] + 1.96 * tot_sd_ood

      eq_models[[j]] <- list(intercept = inter_mpm, coefs = coefs_mpm)
    }
  }

  predict_rhs_fn <- function(X_state_mat) {
    Phi_new <- build_Phi_fn(X_state_mat)
    out_d <- matrix(0.0, nrow = nrow(Phi_new), ncol = m)
    for (j in seq_len(m)) {
      out_d[, j] <- eq_models[[j]]$intercept + as.numeric(Phi_new %*% eq_models[[j]]$coefs)
    }
    out_d
  }

  list(
    method_id = method_id,
    discovered = discovered,
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

