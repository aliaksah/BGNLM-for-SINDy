# ==============================================================================
# R/05_symbolic_canonicalize.R
# Symbolic Canonicalization & Algebraic Equivalence Parser for FBMS, SINDy, & SR
# ==============================================================================

#' Normalize variable names and FBMS / SR / PySINDy function aliases
preprocess_math_string <- function(expr_str) {
  if (is.null(expr_str) || is.na(expr_str)) return("")
  s <- trimws(as.character(expr_str))
  if (!nzchar(s)) return("")

  # Strip outer I(...) wrapper if present
  s <- gsub("^I\\((.*)\\)$", "\\1", s)
  # Normalize x.1 -> x1, x_1 -> x1
  s <- gsub("x\\.([0-9]+)", "x\\1", s)
  s <- gsub("x_([0-9]+)", "x\\1", s)
  # Python power ** to ^
  s <- gsub("\\*\\*", "^", s)
  # Normalize function aliases
  s <- gsub("\\bsin_rad\\b", "sin", s)
  s <- gsub("\\bcos_rad\\b", "cos", s)
  s <- gsub("\\bsin_deg\\b", "sin", s)
  s <- gsub("\\bcos_deg\\b", "cos", s)
  s <- gsub("\\bsqroot\\b", "sqrt", s)
  s <- gsub("\\bp05\\b", "sqrt", s)
  s <- gsub("\\bsquare\\b", "p2", s)
  s <- gsub("\\bcube\\b", "p3", s)

  # PySINDy space-separated monomials e.g. "x1 x2" or "x1^2 x2"
  s <- gsub("([a-zA-Z0-9_\\)\\^])\\s+([a-zA-Z])", "\\1*\\2", s)
  # Remove all whitespace
  s <- gsub("\\s+", "", s)
  s
}

#' Convert an R language AST node into a list of additive monomials,
#' where each monomial is list(coef = numeric, num = character_vec, den = character_vec).
#' Bounded by `depth` and `max_terms` to prevent combinatorial explosion on deeply nested features.
ast_to_monomials <- function(node, depth = 0L, max_terms = 128L) {
  if (is.numeric(node)) {
    return(list(list(coef = as.numeric(node), num = character(0), den = character(0))))
  }
  if (is.symbol(node)) {
    nm <- as.character(node)
    if (nm == "pi") {
      return(list(list(coef = pi, num = character(0), den = character(0))))
    }
    return(list(list(coef = 1.0, num = nm, den = character(0))))
  }
  if (!is.call(node)) {
    return(list(list(coef = 1.0, num = deparse(node), den = character(0))))
  }
  if (depth > 14L) {
    expr_txt <- paste(deparse(node, width.cutoff = 500L), collapse = "")
    return(list(list(coef = 1.0, num = expr_txt, den = character(0))))
  }

  op <- as.character(node[[1]])
  args <- as.list(node)[-1]

  # Parentheses
  if (op == "(") {
    return(ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms))
  }

  # Unary plus / minus
  if (op == "+" && length(args) == 1L) {
    return(ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms))
  }
  if (op == "-" && length(args) == 1L) {
    m1 <- ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms)
    return(lapply(m1, function(m) { m$coef <- -m$coef; m }))
  }

  # Binary addition / subtraction
  if (op == "+" && length(args) == 2L) {
    combined <- c(
      ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms),
      ast_to_monomials(args[[2]], depth = depth + 1L, max_terms = max_terms)
    )
    if (length(combined) > max_terms) {
      expr_txt <- paste(deparse(node, width.cutoff = 500L), collapse = "")
      return(list(list(coef = 1.0, num = expr_txt, den = character(0))))
    }
    return(combined)
  }
  if (op == "-" && length(args) == 2L) {
    m2 <- lapply(ast_to_monomials(args[[2]], depth = depth + 1L, max_terms = max_terms), function(m) { m$coef <- -m$coef; m })
    combined <- c(ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms), m2)
    if (length(combined) > max_terms) {
      expr_txt <- paste(deparse(node, width.cutoff = 500L), collapse = "")
      return(list(list(coef = 1.0, num = expr_txt, den = character(0))))
    }
    return(combined)
  }

  # Multiplication (distributes over sums, guarded against combinatorial blowup)
  if (op == "*" && length(args) == 2L) {
    left_list <- ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms)
    right_list <- ast_to_monomials(args[[2]], depth = depth + 1L, max_terms = max_terms)
    n_prod <- as.numeric(length(left_list)) * as.numeric(length(right_list))
    if (n_prod > max_terms) {
      expr_txt <- paste(deparse(node, width.cutoff = 500L), collapse = "")
      return(list(list(coef = 1.0, num = expr_txt, den = character(0))))
    }
    out <- vector("list", as.integer(n_prod))
    idx_out <- 0L
    for (ml in left_list) {
      for (mr in right_list) {
        idx_out <- idx_out + 1L
        out[[idx_out]] <- list(
          coef = ml$coef * mr$coef,
          num = c(ml$num, mr$num),
          den = c(ml$den, mr$den)
        )
      }
    }
    return(out)
  }

  # Division
  if (op == "/" && length(args) == 2L) {
    left_list <- ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms)
    right_node <- args[[2]]
    # Check if denominator is purely numeric
    if (is.numeric(right_node)) {
      d_val <- as.numeric(right_node)
      return(lapply(left_list, function(ml) { ml$coef <- ml$coef / d_val; ml }))
    }
    den_str <- canonicalize_denominator(right_node)
    return(lapply(left_list, function(ml) {
      ml$den <- c(ml$den, den_str)
      ml
    }))
  }

  # Exponentiation
  if (op == "^" && length(args) == 2L) {
    base_node <- args[[1]]
    exp_node <- args[[2]]
    if (is.numeric(exp_node)) {
      p <- as.numeric(exp_node)
      if (p == 0) {
        return(list(list(coef = 1.0, num = character(0), den = character(0))))
      }
      if (p == 0.5) {
        arg_c <- canonicalize_single_expression(paste(deparse(base_node, width.cutoff = 500L), collapse = ""))
        return(list(list(coef = 1.0, num = paste0("sqrt(", arg_c, ")"), den = character(0))))
      }
      if (p == round(p) && p > 0 && p <= 6) {
        base_list <- ast_to_monomials(base_node, depth = depth + 1L, max_terms = max_terms)
        if ((as.numeric(length(base_list)) ^ p) > max_terms) {
          expr_txt <- paste(deparse(node, width.cutoff = 500L), collapse = "")
          return(list(list(coef = 1.0, num = expr_txt, den = character(0))))
        }
        res <- base_list
        if (p >= 2) {
          for (k in 2:p) {
            n_next <- length(res) * length(base_list)
            next_res <- vector("list", n_next)
            idx_n <- 0L
            for (ml in res) {
              for (mr in base_list) {
                idx_n <- idx_n + 1L
                next_res[[idx_n]] <- list(
                  coef = ml$coef * mr$coef,
                  num = c(ml$num, mr$num),
                  den = c(ml$den, mr$den)
                )
              }
            }
            res <- next_res
          }
        }
        return(res)
      }
    }
  }

  # FBMS custom operator transforms
  if (op == "p2" && length(args) == 1L) {
    return(ast_to_monomials(call("^", args[[1]], 2), depth = depth + 1L, max_terms = max_terms))
  }
  if (op == "p3" && length(args) == 1L) {
    return(ast_to_monomials(call("^", args[[1]], 3), depth = depth + 1L, max_terms = max_terms))
  }
  if (op == "p0" && length(args) == 1L) {
    return(list(list(coef = 1.0, num = character(0), den = character(0))))
  }
  if (op == "pm1" && length(args) == 1L) {
    den_s <- canonicalize_denominator(args[[1]])
    return(list(list(coef = 1.0, num = character(0), den = den_s)))
  }
  if (op == "pm2" && length(args) == 1L) {
    den_s <- canonicalize_denominator(call("^", args[[1]], 2))
    return(list(list(coef = 1.0, num = character(0), den = den_s)))
  }
  if (op == "inv_pos" && length(args) == 1L) {
    arg_s <- canonicalize_single_expression(paste(deparse(args[[1]], width.cutoff = 500L), collapse = ""))
    return(list(list(coef = 1.0, num = character(0), den = paste0("1+", arg_s))))
  }
  if (op == "inv_sq" && length(args) == 1L) {
    sq_m <- ast_to_monomials(call("^", args[[1]], 2), depth = depth + 1L, max_terms = max_terms)[[1]]
    sq_s <- format_monomial_atom(sq_m)
    return(list(list(coef = 1.0, num = character(0), den = paste0("1+", sq_s))))
  }
  if (op == "hill_2" && length(args) == 1L) {
    sq_m <- ast_to_monomials(call("^", args[[1]], 2), depth = depth + 1L, max_terms = max_terms)[[1]]
    sq_s <- format_monomial_atom(sq_m)
    return(list(list(coef = 1.0, num = sq_m$num, den = paste0("1+", sq_s))))
  }
  if (op == "abs" && length(args) == 1L) {
    return(ast_to_monomials(args[[1]], depth = depth + 1L, max_terms = max_terms))
  }

  # General unary functions: sin, cos, exp, sqrt, sigmoid, etc.
  arg_str <- paste(deparse(args[[1]]), collapse = "")
  arg_canon <- canonicalize_single_expression(arg_str)
  fn_token <- paste0(op, "(", arg_canon, ")")
  list(list(coef = 1.0, num = fn_token, den = character(0)))
}

#' Format denominator node into canonical string such as "1+x1" or "1+x3^2"
canonicalize_denominator <- function(node) {
  mons <- ast_to_monomials(node)
  atoms <- character(0)
  has_one <- FALSE
  for (m in mons) {
    a <- format_monomial_atom(m)
    if (a == "1") {
      has_one <- TRUE
    } else if (nzchar(a)) {
      atoms <- c(atoms, a)
    }
  }
  atoms <- sort(unique(atoms))
  if (has_one) {
    return(paste(c("1", atoms), collapse = "+"))
  }
  if (length(atoms) == 0L) return("1")
  paste(atoms, collapse = "+")
}

#' Format a single monomial (num, den) into its canonical basis atom string
format_monomial_atom <- function(m) {
  if (abs(m$coef) < 1e-12) return("")

  # Expand any factor of the form v^k into k copies of v
  expanded_num <- character(0)
  for (f in m$num) {
    if (grepl("^([a-zA-Z0-9_]+)\\^([0-9]+)$", f)) {
      base_v <- sub("^([a-zA-Z0-9_]+)\\^([0-9]+)$", "\\1", f)
      pwr <- as.integer(sub("^([a-zA-Z0-9_]+)\\^([0-9]+)$", "\\2", f))
      expanded_num <- c(expanded_num, rep(base_v, pwr))
    } else if (f != "1" && nzchar(f)) {
      expanded_num <- c(expanded_num, f)
    }
  }

  if (length(expanded_num) == 0L) {
    num_str <- "1"
  } else {
    tbl <- table(sort(expanded_num))
    parts <- character(0)
    for (nm in sort(names(tbl))) {
      cnt <- tbl[[nm]]
      if (cnt == 1L) {
        parts <- c(parts, nm)
      } else {
        parts <- c(parts, paste0(nm, "^", cnt))
      }
    }
    num_str <- paste(parts, collapse = "*")
  }

  den_clean <- m$den[nzchar(m$den) & m$den != "1"]
  if (length(den_clean) == 0L) {
    return(num_str)
  }
  den_sorted <- sort(den_clean)
  den_str <- paste(paste0("(", den_sorted, ")"), collapse = "*")
  paste0(num_str, "/", den_str)
}

#' Canonicalize a single feature / basis expression string
canonicalize_single_expression <- function(expr_str) {
  atoms <- extract_canonical_atoms(expr_str)
  if (length(atoms) == 0L) return("1")
  paste(atoms, collapse = "+")
}

#' Extract canonical basis atoms (excluding intercept "1") from a feature or equation string
extract_canonical_atoms <- function(expr_strings, include_intercept = FALSE) {
  if (length(expr_strings) == 0L) return(character(0))
  all_atoms <- character(0)

  for (raw_s in expr_strings) {
    s <- preprocess_math_string(raw_s)
    if (!nzchar(s) || s %in% c("1", "(Intercept)", "Intercept")) {
      if (include_intercept) all_atoms <- c(all_atoms, "1")
      next
    }
    if (nchar(s) > 500L) {
      all_atoms <- c(all_atoms, substr(s, 1L, 500L))
      next
    }
    parsed <- tryCatch(parse(text = s)[[1]], error = function(e) NULL)
    if (is.null(parsed)) {
      all_atoms <- c(all_atoms, s)
      next
    }
    mons <- tryCatch(ast_to_monomials(parsed), error = function(e) NULL)
    if (is.null(mons)) {
      all_atoms <- c(all_atoms, s)
      next
    }
    for (m in mons) {
      atom <- format_monomial_atom(m)
      if (nzchar(atom)) {
        if (atom == "1" && !include_intercept) next
        all_atoms <- c(all_atoms, atom)
      }
    }
  }
  sort(unique(all_atoms))
}

#' Prune negligible basis terms based on relative L2 energy ||beta_k g_k(X)||_2 / ||Y||_2 < tol_rel
prune_negligible_terms <- function(atom_strings, coefs, X, y, tol_rel = 1e-3) {
  if (length(atom_strings) == 0L) return(character(0))
  y_norm <- sqrt(sum(y^2))
  if (!is.finite(y_norm) || y_norm < 1e-12) y_norm <- 1.0

  keep <- logical(length(atom_strings))
  for (k in seq_along(atom_strings)) {
    atom <- atom_strings[k]
    if (atom %in% c("1", "(Intercept)", "Intercept")) {
      keep[k] <- FALSE
      next
    }
    beta_k <- if (length(coefs) >= k && is.finite(coefs[k])) coefs[k] else 1.0
    if (abs(beta_k) < 1e-10) {
      keep[k] <- FALSE
      next
    }
    gk <- evaluate_basis_atom(atom, X)
    term_norm <- abs(beta_k) * sqrt(sum(gk^2))
    keep[k] <- (term_norm / y_norm) >= tol_rel
  }
  extract_canonical_atoms(atom_strings[keep])
}

#' Count total mathematical operations and variables (Operation Count OC) in discovered atoms
compute_operation_count <- function(atom_strings) {
  atoms <- extract_canonical_atoms(atom_strings)
  if (length(atoms) == 0L) return(0L)
  total_oc <- 0L
  for (a in atoms) {
    # Count variables x1..xm
    vars <- gregexpr("x[0-9]+", a)[[1]]
    n_vars <- if (vars[1] == -1L) 0L else length(vars)
    # Count binary operators *, /, +, ^
    ops <- gregexpr("[\\*\\/\\+\\^]", a)[[1]]
    n_ops <- if (ops[1] == -1L) 0L else length(ops)
    # Count unary functions sin, cos, exp, sqrt, sigmoid
    fns <- gregexpr("(sin|cos|exp|sqrt|sigmoid)", a)[[1]]
    n_fns <- if (fns[1] == -1L) 0L else length(fns)
    total_oc <- total_oc + n_vars + n_ops + n_fns
  }
  as.integer(total_oc)
}

