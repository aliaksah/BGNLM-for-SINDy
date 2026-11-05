# ==============================================================================
# R/06_metrics_consistency.R
# Selection Consistency, Stability (Jaccard & Nogueira), PIP Calibration, & Rollout
# ==============================================================================

if (!exists("extract_canonical_atoms", mode = "function")) {
  if (file.exists("R/05_symbolic_canonicalize.R")) {
    source("R/05_symbolic_canonicalize.R")
  } else if (file.exists("../../R/05_symbolic_canonicalize.R")) {
    source("../../R/05_symbolic_canonicalize.R")
  }
}

#' Compute structural identification metrics for a single equation or full ODE system
#' @param discovered_list List of character vectors of discovered basis strings per equation
#' @param ground_truth_list List of character vectors of true basis strings per equation
compute_structure_metrics <- function(discovered_list, ground_truth_list) {
  m <- length(ground_truth_list)
  stopifnot(length(discovered_list) == m)

  # tag each canonical atom with its equation index "eqj::atom" for exact system-wide matching
  S_hat <- character(0)
  S_star <- character(0)
  tpr_eq <- numeric(m)
  fdr_eq <- numeric(m)
  exact_eq <- logical(m)
  oc_total <- 0L

  for (j in seq_len(m)) {
    hat_j <- extract_canonical_atoms(discovered_list[[j]])
    star_j <- extract_canonical_atoms(ground_truth_list[[j]])
    oc_total <- oc_total + compute_operation_count(hat_j)

    tp_j <- length(intersect(hat_j, star_j))
    fp_j <- length(setdiff(hat_j, star_j))
    fn_j <- length(setdiff(star_j, hat_j))

    tpr_eq[j] <- if (length(star_j) > 0L) tp_j / length(star_j) else 1.0
    fdr_eq[j] <- if (length(hat_j) > 0L) fp_j / length(hat_j) else 0.0
    exact_eq[j] <- (fp_j == 0L && fn_j == 0L)

    if (length(hat_j) > 0L) S_hat <- c(S_hat, paste0("eq", j, "::", hat_j))
    if (length(star_j) > 0L) S_star <- c(S_star, paste0("eq", j, "::", star_j))
  }

  tp <- length(intersect(S_hat, S_star))
  fp <- length(setdiff(S_hat, S_star))
  fn <- length(setdiff(S_star, S_hat))

  tpr <- if (length(S_star) > 0L) tp / length(S_star) else 1.0
  fdr <- if (length(S_hat) > 0L) fp / length(S_hat) else 0.0
  precision <- 1.0 - fdr
  f1 <- if ((tpr + precision) > 0) (2.0 * tpr * precision) / (tpr + precision) else 0.0
  shd <- as.integer(fp + fn)
  exact_recovery <- as.integer(fp == 0L && fn == 0L)

  list(
    TPR = tpr,
    FDR = fdr,
    F1 = f1,
    Exact_Recovery = exact_recovery,
    SHD = shd,
    OC = oc_total,
    n_active_terms = length(S_hat),
    TPR_eq_mean = mean(tpr_eq),
    FDR_eq_mean = mean(fdr_eq),
    S_hat = S_hat,
    S_star = S_star
  )
}

#' Compute Posterior Inclusion Probability (PIP) Brier Score across candidate features
#' @param pip_list List of named numeric vectors (names = feature strings, values = PIP in [0,1]) per equation
#' @param ground_truth_list List of true basis atom vectors per equation
compute_pip_brier_score <- function(pip_list, ground_truth_list) {
  if (is.null(pip_list) || length(pip_list) == 0L) return(NA_real_)
  sq_errors <- numeric(0)

  for (j in seq_along(ground_truth_list)) {
    if (j > length(pip_list)) next
    pips <- pip_list[[j]]
    if (is.null(pips) || length(pips) == 0L) next
    star_j <- extract_canonical_atoms(ground_truth_list[[j]])

    # Map candidate feature names to canonical atoms and aggregate max PIP per canonical atom
    canon_names <- vapply(names(pips), function(nm) {
      ca <- extract_canonical_atoms(nm)
      if (length(ca) == 0L) "" else paste(ca, collapse = "+")
    }, character(1))
    valid <- nzchar(canon_names)
    pips <- pips[valid]
    canon_names <- canon_names[valid]

    all_atoms <- union(unique(canon_names), star_j)
    for (a in all_atoms) {
      prob_a <- if (a %in% canon_names) max(pips[canon_names == a], na.rm = TRUE) else 0.0
      prob_a <- min(max(prob_a, 0.0), 1.0)
      truth_a <- as.numeric(a %in% star_j)
      sq_errors <- c(sq_errors, (prob_a - truth_a)^2)
    }
  }
  if (length(sq_errors) == 0L) return(NA_real_)
  mean(sq_errors)
}

#' Compute R^2, clipped R^2, and NRMSE between predicted and true derivatives/trajectories
compute_prediction_metrics <- function(Y_pred, Y_true, lower_95 = NULL, upper_95 = NULL) {
  Y_pred <- as.matrix(Y_pred)
  Y_true <- as.matrix(Y_true)
  m <- ncol(Y_true)

  r2_raw_vec <- numeric(m)
  nrmse_vec <- numeric(m)
  rmse_vec <- numeric(m)

  for (j in seq_len(m)) {
    yp <- Y_pred[, j]
    yt <- Y_true[, j]
    ok <- is.finite(yp) & is.finite(yt)
    if (sum(ok) < 3L) {
      r2_raw_vec[j] <- 0.0
      nrmse_vec[j] <- 1.0
      rmse_vec[j] <- NA_real_
      next
    }
    yp <- yp[ok]; yt <- yt[ok]
    ss_res <- sum((yt - yp)^2)
    ss_tot <- sum((yt - mean(yt))^2)
    rmse_j <- sqrt(mean((yt - yp)^2))
    sd_y <- sd(yt)
    if (!is.finite(sd_y) || sd_y < 1e-12) sd_y <- 1.0

    r2_sse <- 1.0 - (ss_res / max(ss_tot, 1e-12))
    r2_cor <- suppressWarnings(cor(yp, yt)^2)
    if (!is.finite(r2_cor)) r2_cor <- 0.0
    # Store standard R2 (with correlation R2 fallback if constant offset)
    r2_raw_vec[j] <- r2_sse
    nrmse_vec[j] <- rmse_j / sd_y
    rmse_vec[j] <- rmse_j
  }

  r2_raw <- mean(r2_raw_vec, na.rm = TRUE)
  r2_clipped <- mean(pmax(0.0, pmin(1.0, r2_raw_vec)), na.rm = TRUE)
  nrmse <- mean(pmin(nrmse_vec, 10.0), na.rm = TRUE)
  log10_rmse <- log10(max(mean(rmse_vec, na.rm = TRUE), 1e-12))

  # 95% Credible / Prediction Interval Coverage and Winkler Score
  coverage_95 <- NA_real_
  winkler_95 <- NA_real_
  if (!is.null(lower_95) && !is.null(upper_95)) {
    L <- as.matrix(lower_95)
    U <- as.matrix(upper_95)
    if ( all(dim(L) == dim(Y_true)) && all(dim(U) == dim(Y_true)) ) {
      alpha <- 0.05
      inside <- (Y_true >= L) & (Y_true <= U)
      coverage_95 <- mean(inside, na.rm = TRUE)
      width <- pmax(0.0, U - L)
      penalty_low <- (2.0 / alpha) * pmax(0.0, L - Y_true)
      penalty_high <- (2.0 / alpha) * pmax(0.0, Y_true - U)
      sd_scale <- mean(apply(Y_true, 2, sd), na.rm = TRUE)
      if (!is.finite(sd_scale) || sd_scale < 1e-8) sd_scale <- 1.0
      winkler_95 <- mean(width + penalty_low + penalty_high, na.rm = TRUE) / sd_scale
    }
  }

  list(
    R2_raw = r2_raw,
    R2_clipped = r2_clipped,
    NRMSE = nrmse,
    log10_RMSE = log10_rmse,
    Coverage_95 = coverage_95,
    Winkler_95 = winkler_95
  )
}

#' Forward ODE Trajectory Rollout Evaluator (R2_traj, NRMSE_traj, and Valid Prediction Time T_valid)
#' @param predict_rhs_fn Function taking state matrix X (n x m) and returning derivative matrix dX (n x m)
#' @param rollout_spec List with `times`, `X_true`, `x0`, and optional `X_inputs`
evaluate_trajectory_rollout <- function(predict_rhs_fn, rollout_spec, tol_valid = 0.20) {
  times <- rollout_spec$times
  X_true <- as.matrix(rollout_spec$X_true)
  x0 <- as.numeric(rollout_spec$x0)
  m <- ncol(X_true)
  names(x0) <- colnames(X_true)
  X_inputs <- rollout_spec$X_inputs

  bound_max <- 50.0 * max(abs(X_true), 1.0)

  ode_func <- function(t, state, parms) {
    if (any(!is.finite(state)) || max(abs(state)) > bound_max) {
      return(list(rep(0.0, m)))
    }
    state_mat <- matrix(state, nrow = 1L, ncol = m)
    colnames(state_mat) <- colnames(X_true)
    if (!is.null(X_inputs)) {
      # Interpolate exogenous inputs at time t
      idx_t <- which.min(abs(times - t))
      inp_mat <- matrix(X_inputs[idx_t, ], nrow = 1L)
      colnames(inp_mat) <- colnames(X_inputs)
      state_mat <- cbind(state_mat, inp_mat)
    }
    d_state <- tryCatch(
      as.numeric(predict_rhs_fn(state_mat)),
      error = function(e) rep(0.0, m)
    )
    d_state[!is.finite(d_state)] <- 0.0
    d_state <- pmax(pmin(d_state, bound_max), -bound_max)
    list(d_state)
  }

  sol <- tryCatch(
    suppressWarnings(deSolve::ode(
      y = x0,
      times = times,
      func = ode_func,
      parms = NULL,
      method = "rk4"
    )),
    error = function(e) NULL
  )

  if (is.null(sol) || nrow(sol) != length(times)) {
    return(list(R2_traj = 0.0, NRMSE_traj = 1.0, T_valid = 0.0, X_sim = matrix(NA_real_, nrow(X_true), m)))
  }

  X_sim <- as.matrix(sol[, -1, drop = FALSE])
  pred_m <- compute_prediction_metrics(X_sim, X_true)

  # Valid prediction time T_valid until normalized step error exceeds tol_valid
  scale_sd <- apply(X_true, 2, sd)
  scale_sd[scale_sd < 1e-8] <- 1.0
  rel_err_t <- sqrt(rowMeans(sweep(X_sim - X_true, 2, scale_sd, "/")^2))
  exceed_idx <- which(!is.finite(rel_err_t) | rel_err_t > tol_valid)
  if (length(exceed_idx) == 0L) {
    t_valid <- max(times) - min(times)
  } else {
    t_valid <- times[max(1L, exceed_idx[1] - 1L)] - min(times)
  }
  t_horizon <- max(max(times) - min(times), 1e-8)

  list(
    R2_traj = pred_m$R2_clipped,
    R2_traj_raw = pred_m$R2_raw,
    NRMSE_traj = pred_m$NRMSE,
    T_valid = t_valid,
    T_valid_frac = t_valid / t_horizon,
    X_sim = X_sim
  )
}

#' Pairwise Jaccard Stability across R replications
#' @param selected_sets_list List of length R, each a character vector of selected canonical atoms
compute_jaccard_stability <- function(selected_sets_list) {
  R <- length(selected_sets_list)
  if (R < 2L) return(NA_real_)
  j_vals <- numeric(0)
  for (i in 1:(R - 1L)) {
    Si <- unique(selected_sets_list[[i]])
    for (j in (i + 1L):R) {
      Sj <- unique(selected_sets_list[[j]])
      u_len <- length(union(Si, Sj))
      j_ij <- if (u_len == 0L) 1.0 else length(intersect(Si, Sj)) / u_len
      j_vals <- c(j_vals, j_ij)
    }
  }
  mean(j_vals)
}

#' Nogueira's Stability Estimator (Nogueira et al., JMLR 2018) with 95% Confidence Intervals
#' @param selected_sets_list List of length R, each a character vector of selected canonical atoms
#' @param d_total Optional total universe size of candidate features (inferred from union + padding if NULL)
compute_nogueira_stability <- function(selected_sets_list, d_total = NULL) {
  M <- length(selected_sets_list)
  if (M < 2L) {
    return(list(phi = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_))
  }
  universe <- sort(unique(unlist(selected_sets_list)))
  if (length(universe) == 0L) {
    return(list(phi = 1.0, ci_lower = 1.0, ci_upper = 1.0))
  }
  d <- if (is.null(d_total)) max(length(universe) + 10L, 20L) else max(as.integer(d_total), length(universe))

  # Construct binary selection matrix Z (M x d)
  Z <- matrix(0.0, nrow = M, ncol = d)
  for (m in seq_len(M)) {
    idx <- match(unique(selected_sets_list[[m]]), universe)
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0L) Z[m, idx] <- 1.0
  }

  p_hat <- colMeans(Z)
  s_f2 <- (M / (M - 1.0)) * p_hat * (1.0 - p_hat)
  k_bar <- sum(rowSums(Z)) / M
  denom <- (k_bar / d) * (1.0 - k_bar / d)

  if (denom < 1e-12) {
    return(list(phi = 1.0, ci_lower = 1.0, ci_upper = 1.0))
  }
  phi_hat <- 1.0 - (mean(s_f2) / denom)
  phi_hat <- max(min(phi_hat, 1.0), -1.0)

  # Jackknife / influence-function variance estimate (Nogueira et al. 2018 Theorem 5)
  phi_i <- numeric(M)
  for (i in seq_len(M)) {
    phi_i[i] <- (1.0 / denom) * (mean(Z[i, ] * p_hat) - (k_bar * sum(rowSums(Z)[i])) / (d^2) +
                                   (phi_hat / 2.0) * ((2.0 * k_bar * sum(Z[i, ])) / (d^2) - sum(Z[i, ]) / d - k_bar / d + 1.0))
  }
  var_phi <- (4.0 / (M^2)) * sum((phi_i - mean(phi_i))^2)
  se_phi <- sqrt(max(var_phi, 0.0))

  list(
    phi = phi_hat,
    se = se_phi,
    ci_lower = max(-1.0, phi_hat - 1.96 * se_phi),
    ci_upper = min(1.0, phi_hat + 1.96 * se_phi)
  )
}

