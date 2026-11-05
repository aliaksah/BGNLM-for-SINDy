# ==============================================================================
# R/04b_sr_baselines.R
# Symbolic Regression Baselines (SR_PySR, SR_Operon) via Python Bridge
# ==============================================================================

if (!exists("extract_canonical_atoms", mode = "function")) {
  if (file.exists("R/05_symbolic_canonicalize.R")) {
    source("R/05_symbolic_canonicalize.R")
  } else if (file.exists("../../R/05_symbolic_canonicalize.R")) source("../../R/05_symbolic_canonicalize.R")
}

#' Locate Python interpreter in .venv or system PATH (preserving venv symlink)
find_python_bin <- function() {
  candidates <- c(
    ".venv/bin/python",
    "../../.venv/bin/python",
    Sys.which("python3"),
    Sys.which("python")
  )
  for (cand in candidates) {
    if (nzchar(cand) && file.exists(cand)) {
      return(file.path(normalizePath(dirname(cand), mustWork = TRUE), basename(cand)))
    }
  }
  stop("No Python interpreter found for Symbolic Regression bridge.")
}

#' Locate python/run_sr_pysindy.py script
find_sr_script <- function() {
  candidates <- c("python/run_sr_pysindy.py", "../../python/run_sr_pysindy.py")
  for (cand in candidates) {
    if (file.exists(cand)) {
      return(normalizePath(cand))
    }
  }
  stop("Cannot find python/run_sr_pysindy.py")
}

#' Run Symbolic Regression Baseline (SR_PySR or SR_Operon)
run_sr_baseline <- function(
    method_id = c("SR_PySR", "SR_Operon"),
    task_data) {
  method_id <- match.arg(method_id)
  X_tr <- as.matrix(task_data$X_train)
  Y_tr <- as.matrix(task_data$Y_train)
  X_in <- as.matrix(task_data$X_test_in)
  X_ood <- as.matrix(task_data$X_test_ood)
  m <- ncol(Y_tr)
  col_names <- colnames(X_tr)
  if (is.null(col_names)) col_names <- paste0("x", seq_len(ncol(X_tr)))

  py_bin <- find_python_bin()
  sr_script <- find_sr_script()

  tmp_in <- tempfile(pattern = "sr_in_", fileext = ".json")
  tmp_out <- tempfile(pattern = "sr_out_", fileext = ".json")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

  # Write input JSON using jsonlite (installed with R packages)
  payload <- list(
    method_id = method_id,
    seed = as.integer(if (!is.null(task_data$seed)) task_data$seed else 42L),
    col_names = as.list(col_names),
    X_train = X_tr,
    Y_train = Y_tr,
    X_test_in = X_in,
    X_test_ood = X_ood
  )
  jsonlite::write_json(payload, tmp_in, auto_unbox = TRUE, digits = 10)

  cmd_status <- system2(
    py_bin,
    args = c(shQuote(sr_script), "--input", shQuote(tmp_in), "--output", shQuote(tmp_out)),
    stdout = FALSE,
    stderr = FALSE
  )
  if (cmd_status != 0L || !file.exists(tmp_out)) {
    stop("Python SR bridge failed for method: ", method_id)
  }

  res_py <- jsonlite::read_json(tmp_out, simplifyVector = TRUE)

  discovered <- vector("list", m)
  pip_list <- vector("list", m)
  eq_models <- vector("list", m)

  for (j in seq_len(m)) {
    raw_atoms <- if (length(res_py$discovered) >= j) as.character(res_py$discovered[[j]]) else character(0)
    canon_j <- extract_canonical_atoms(raw_atoms)
    discovered[[j]] <- canon_j
    # Symbolic regression yields point-estimate Pareto trees without calibrated PIPs
    pip_list[[j]] <- NULL

    inter_j <- as.numeric(res_py$intercepts[j])
    coef_list_j <- res_py$coefs[[j]]
    if (is.data.frame(res_py$coefs)) {
      coef_vec <- as.numeric(res_py$coefs[j, ])
      names(coef_vec) <- colnames(res_py$coefs)
      coef_vec <- coef_vec[!is.na(coef_vec) & abs(coef_vec) > 0]
    } else if (is.list(coef_list_j) && length(coef_list_j) > 0L) {
      coef_vec <- unlist(coef_list_j)
    } else {
      coef_vec <- numeric(0)
    }
    eq_models[[j]] <- list(intercept = inter_j, coefs = coef_vec)
  }

  pred_train <- as.matrix(res_py$pred_train)
  pred_in <- as.matrix(res_py$pred_in)
  pred_ood <- as.matrix(res_py$pred_ood)

  # Symbolic regression does not provide Bayesian credible intervals
  na_in <- matrix(NA_real_, nrow = nrow(X_in), ncol = m)
  na_ood <- matrix(NA_real_, nrow = nrow(X_ood), ncol = m)

  predict_rhs_fn <- function(X_state_mat) {
    X_df_new <- as.data.frame(X_state_mat)
    colnames(X_df_new) <- col_names
    out_d <- matrix(0.0, nrow = nrow(X_df_new), ncol = m)
    for (j in seq_len(m)) {
      mod_j <- eq_models[[j]]
      val <- rep(mod_j$intercept, nrow(X_df_new))
      if (length(mod_j$coefs) > 0L) {
        for (nm in names(mod_j$coefs)) {
          val <- val + mod_j$coefs[[nm]] * evaluate_basis_atom(nm, X_df_new)
        }
      }
      out_d[, j] <- val
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
    lower_in = na_in,
    upper_in = na_in,
    lower_ood = na_ood,
    upper_ood = na_ood,
    predict_rhs_fn = predict_rhs_fn
  )
}
