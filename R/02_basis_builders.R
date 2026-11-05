# ==============================================================================
# R/02_basis_builders.R
# Predefined Design Matrix Builders: Poly(d=3), Oracle, and Overcomplete
# ==============================================================================

#' Generate all polynomial monomials up to degree `max_deg` for state matrix `X`
build_poly_library <- function(X, max_deg = 3L) {
  X <- as.matrix(X)
  col_names <- colnames(X)
  if (is.null(col_names)) {
    col_names <- paste0("x", seq_len(ncol(X)))
    colnames(X) <- col_names
  }
  m <- ncol(X)
  n <- nrow(X)

  feats_list <- list()
  feat_names <- character(0)

  # Degree 1
  for (i in seq_len(m)) {
    feats_list[[length(feats_list) + 1L]] <- X[, i]
    feat_names <- c(feat_names, col_names[i])
  }

  # Degree 2
  if (max_deg >= 2L) {
    for (i in seq_len(m)) {
      for (j in i:m) {
        feats_list[[length(feats_list) + 1L]] <- X[, i] * X[, j]
        if (i == j) {
          feat_names <- c(feat_names, paste0(col_names[i], "^2"))
        } else {
          feat_names <- c(feat_names, paste0(col_names[i], "*", col_names[j]))
        }
      }
    }
  }

  # Degree 3
  if (max_deg >= 3L) {
    for (i in seq_len(m)) {
      for (j in i:m) {
        for (k in j:m) {
          feats_list[[length(feats_list) + 1L]] <- X[, i] * X[, j] * X[, k]
          idx <- c(i, j, k)
          tbl <- table(idx)
          parts <- character(0)
          for (v_idx in as.integer(names(tbl))) {
            pwr <- tbl[[as.character(v_idx)]]
            if (pwr == 1L) {
              parts <- c(parts, col_names[v_idx])
            } else {
              parts <- c(parts, paste0(col_names[v_idx], "^", pwr))
            }
          }
          feat_names <- c(feat_names, paste(parts, collapse = "*"))
        }
      }
    }
  }

  Phi <- do.call(cbind, feats_list)
  colnames(Phi) <- feat_names
  Phi
}

#' Evaluate a named symbolic basis atom on state matrix `X`
evaluate_basis_atom <- function(atom_str, X) {
  df <- as.data.frame(X)
  # Safe helpers for fractional power and rationals on states
  env <- list2env(as.list(df), parent = baseenv())
  env$sqrt <- function(z) base::sqrt(abs(z))
  # Replace "1/(1+x" or "+x" denominators with safe abs() if negative states occur
  val <- tryCatch(
    eval(parse(text = atom_str), envir = env),
    error = function(e) rep(0.0, nrow(X))
  )
  val[!is.finite(val)] <- 0.0
  as.numeric(val)
}

#' Build Known Correct (Oracle) Basis: Degree-3 Polynomials + Exact Ground-Truth Nonlinear Terms
build_oracle_library <- function(X, extra_atoms = character(0), max_deg = 3L) {
  Phi_poly <- build_poly_library(X, max_deg = max_deg)
  existing_names <- colnames(Phi_poly)
  extra_atoms <- setdiff(unique(extra_atoms), existing_names)
  if (length(extra_atoms) == 0L) {
    return(Phi_poly)
  }
  extra_cols <- lapply(extra_atoms, function(a) evaluate_basis_atom(a, X))
  Phi_extra <- do.call(cbind, extra_cols)
  colnames(Phi_extra) <- extra_atoms
  cbind(Phi_poly, Phi_extra)
}

#' Build Brute-Force Overcomplete Library (q > 60 columns: Poly(3) + Trig + Exp + Rational + Fractional + Cross-Terms)
build_overcomplete_library <- function(X, extra_atoms = character(0)) {
  X <- as.matrix(X)
  col_names <- colnames(X)
  if (is.null(col_names)) {
    col_names <- paste0("x", seq_len(ncol(X)))
    colnames(X) <- col_names
  }
  m <- ncol(X)
  Phi_poly <- build_poly_library(X, max_deg = 3L)

  extra_list <- list()
  extra_names <- character(0)
  add_col <- function(vec, nm) {
    if (!nm %in% c(colnames(Phi_poly), extra_names)) {
      vec[!is.finite(vec)] <- 0.0
      extra_list[[length(extra_list) + 1L]] <<- as.numeric(vec)
      extra_names <<- c(extra_names, nm)
    }
  }

  # Unary non-polynomials for every state variable
  for (i in seq_len(m)) {
    xi <- X[, i]
    nm <- col_names[i]
    add_col(sin(xi), paste0("sin(", nm, ")"))
    add_col(cos(xi), paste0("cos(", nm, ")"))
    add_col(sin(2 * xi), paste0("sin(2*", nm, ")"))
    add_col(cos(2 * xi), paste0("cos(2*", nm, ")"))
    add_col(sin(3 * xi), paste0("sin(3*", nm, ")"))
    add_col(cos(3 * xi), paste0("cos(3*", nm, ")"))
    add_col(exp(-abs(xi)), paste0("exp(-abs(", nm, "))"))
    add_col(1.0 / (1.0 + abs(xi)), paste0("1/(1+", nm, ")"))
    add_col(1.0 / (1.0 + xi^2), paste0("1/(1+", nm, "^2)"))
    add_col(xi / (1.0 + abs(xi)), paste0(nm, "/(1+", nm, ")"))
    add_col((xi^2) / (1.0 + xi^2), paste0(nm, "^2/(1+", nm, "^2)"))
    add_col(sqrt(abs(xi)), paste0("sqrt(", nm, ")"))
    add_col(1.0 / (1.0 + exp(-pmin(pmax(xi, -20), 20))), paste0("sigmoid(", nm, ")"))
    add_col(xi^4, paste0(nm, "^4"))
  }

  # Cross-interactions between states and non-polynomials to ensure > 60 collinear features
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      xi <- X[, i]; xj <- X[, j]
      ni <- col_names[i]; nj <- col_names[j]
      add_col(xi * sin(xj), paste0(ni, "*sin(", nj, ")"))
      add_col(xi * cos(xj), paste0(ni, "*cos(", nj, ")"))
      add_col((xi * xj) / (1.0 + abs(xi)), paste0(ni, "*", nj, "/(1+", ni, ")"))
      add_col((xi * xj^2) / (1.0 + xj^2), paste0(ni, "*", nj, "^2/(1+", nj, "^2)"))
      add_col(xi * exp(-abs(xj)), paste0(ni, "*exp(-abs(", nj, "))"))
      add_col(sin(xi) * cos(xj), paste0("sin(", ni, ")*cos(", nj, ")"))
    }
  }

  for (a in extra_atoms) {
    add_col(evaluate_basis_atom(a, X), a)
  }

  Phi_extra <- do.call(cbind, extra_list)
  colnames(Phi_extra) <- extra_names
  cbind(Phi_poly, Phi_extra)
}

