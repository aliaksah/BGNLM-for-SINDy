#!/usr/bin/env Rscript
# ==============================================================================
# scripts/install_dependencies.R
# Installs Required R Packages into Project-Local .r_lib/ & Configures .venv/
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
use_global <- "--global" %in% args

# Resolve repo root
repo_root <- normalizePath(".", mustWork = TRUE)
local_r_lib <- file.path(repo_root, ".r_lib")

if (!use_global) {
  if (!dir.exists(local_r_lib)) {
    dir.create(local_r_lib, recursive = TRUE, showWarnings = FALSE)
  }
  local_r_lib <- normalizePath(local_r_lib, mustWork = TRUE)
  target_lib <- local_r_lib
  .libPaths(unique(c(local_r_lib, .Library, .libPaths())))
  Sys.setenv(R_LIBS_USER = local_r_lib, R_LIBS = local_r_lib)
  message("[ENV] Installing R packages into project-local library: ", target_lib)
} else {
  user_lib <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(user_lib)) {
    user_lib <- path.expand(file.path(
      "~", "R",
      paste0(R.version$platform, "-library"),
      paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1])
    ))
  }
  if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  target_lib <- user_lib
  .libPaths(unique(c(target_lib, .libPaths())))
  message("[ENV] Installing R packages into user library: ", target_lib)
}

n_cores <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 2L)
if (is.na(n_cores) || n_cores < 2L) n_cores <- 1L
n_cpus <- min(4L, max(1L, n_cores - 1L))
options(Ncpus = n_cpus)
message("[ENV] Using Ncpus = ", n_cpus, " for parallel source package compilation.")

# Helper: metadata-only check (NEVER call loadNamespace() before install.packages!)
installed_names <- function(lib_loc) {
  rownames(installed.packages(lib.loc = c(lib_loc, .Library)))
}

# Ensure modern rlang / tidyverse core dependencies exist in target_lib so they shadow older system libs
core_infra <- c("rlang", "cli", "glue", "lifecycle", "vctrs", "remotes")
have_in_target <- rownames(installed.packages(lib.loc = target_lib))
need_infra <- setdiff(core_infra, have_in_target)
if (length(need_infra) > 0L) {
  message("[INSTALL] Bootstrapping core R infrastructure in ", target_lib, ": ",
          paste(need_infra, collapse = ", "))
  install.packages(
    need_infra,
    lib = target_lib,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = n_cpus
  )
}

# CRAN packages needed
required_cran_pkgs <- c(
  "remotes", "deSolve", "glmnet", "signal", "data.table", "yaml", "jsonlite",
  "lme4", "lmerTest", "glmmTMB", "emmeans", "broom.mixed",
  "ggplot2", "dplyr", "tidyr", "reticulate", "R.utils", "testthat"
)

have_all <- installed_names(target_lib)
missing_cran <- setdiff(required_cran_pkgs, have_all)

if (length(missing_cran) > 0L) {
  message("[INSTALL] Installing missing CRAN packages (", length(missing_cran), "): ",
          paste(missing_cran, collapse = ", "))
  install.packages(
    missing_cran,
    lib = target_lib,
    repos = "https://cloud.r-project.org",
    dependencies = c("Depends", "Imports", "LinkingTo"),
    Ncpus = n_cpus
  )
}

# Install / update FBMS specifically from GitHub branch: jonlachmann/FBMS@JCS
message("[INSTALL] Installing FBMS from GitHub: jonlachmann/FBMS@JCS into ", target_lib, "...")
fbms_ok <- tryCatch({
  remotes::install_github(
    "jonlachmann/FBMS@JCS",
    lib = target_lib,
    dependencies = c("Depends", "Imports", "LinkingTo"),
    upgrade = "never",
    force = TRUE,
    quiet = FALSE
  )
  "FBMS" %in% rownames(installed.packages(lib.loc = target_lib))
}, error = function(e) {
  message("[WARN] Standard FBMS@JCS install failed: ", e$message)
  FALSE
})

# Fallback if tolerance (which imports plotly) fails on missing RHEL8 curl/openssl headers
if (!fbms_ok) {
  message("[INSTALL] Attempting fallback: installing FBMS runtime deps + tolerance shim if needed...")
  fbms_deps <- c("Rcpp", "RcppArmadillo", "fastglm", "GenSA", "r2r", "BAS")
  need_fbms_deps <- setdiff(fbms_deps, installed_names(target_lib))
  if (length(need_fbms_deps) > 0L) {
    install.packages(
      need_fbms_deps,
      lib = target_lib,
      repos = "https://cloud.r-project.org",
      dependencies = c("Depends", "Imports", "LinkingTo"),
      Ncpus = n_cpus
    )
  }
  if (!("tolerance" %in% installed_names(target_lib))) {
    tol_src <- file.path(tempdir(), "tolerance_shim")
    dir.create(file.path(tol_src, "R"), recursive = TRUE, showWarnings = FALSE)
    writeLines(c(
      "Package: tolerance",
      "Version: 2.0.0",
      "Title: Minimal tolerance shim for FBMS",
      "Author: Shim",
      "Maintainer: Shim <noreply@example.com>",
      "Description: Minimal F1 stub when plotly system headers are unavailable.",
      "License: GPL-2"
    ), file.path(tol_src, "DESCRIPTION"))
    writeLines("export(F1)", file.path(tol_src, "NAMESPACE"))
    writeLines("F1 <- function(...) stop('tolerance::F1 unavailable in minimal shim')", file.path(tol_src, "R", "F1.R"))
    install.packages(tol_src, repos = NULL, type = "source", lib = target_lib)
  }
  remotes::install_github(
    "jonlachmann/FBMS@JCS",
    lib = target_lib,
    dependencies = FALSE,
    upgrade = "never",
    force = TRUE,
    quiet = FALSE
  )
}

required_r_pkgs <- c("FBMS", required_cran_pkgs)

# Re-verify every required package loads cleanly
still_broken <- character(0)
for (pkg in required_r_pkgs) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, lib.loc = c(target_lib, .Library)))
    TRUE
  }, error = function(e) {
    message("[ERROR] Failed loading ", pkg, ": ", e$message)
    FALSE
  })
  if (!ok) still_broken <- c(still_broken, pkg)
}

if (length(still_broken) > 0L) {
  stop("Failed to install/load required R packages: ", paste(still_broken, collapse = ", "))
}

message("[OK] All ", length(required_r_pkgs), " required R packages installed and verified in ",
        target_lib, " (FBMS version: ", as.character(packageVersion("FBMS")), ").")

# Set up project-local Python .venv if missing or incomplete
venv_py <- file.path(repo_root, ".venv", "bin", "python")
py_core <- c("numpy", "scipy", "scikit-learn", "sympy", "pysindy", "pandas", "pyyaml")
py_optional <- c("pysr", "pyoperon")

if (!file.exists(venv_py)) {
  uv_bin <- Sys.which("uv")
  if (!nzchar(uv_bin) && file.exists(path.expand("~/.local/bin/uv"))) {
    uv_bin <- path.expand("~/.local/bin/uv")
  }
  if (nzchar(uv_bin)) {
    message("[INSTALL] Creating .venv via uv...")
    system2(uv_bin, c("venv", ".venv", "--python", "3.11"))
    system2(uv_bin, c("pip", "install", "--python", venv_py, py_core))
    try(system2(uv_bin, c("pip", "install", "--python", venv_py, py_optional)), silent = TRUE)
  } else {
    py3 <- Sys.which("python3")
    if (!nzchar(py3)) stop("No python3 interpreter found on PATH to create .venv.")
    message("[INSTALL] Creating .venv via python3 -m venv...")
    system2(py3, c("-m", "venv", ".venv"))
    system2(".venv/bin/pip", c("install", "--upgrade", "pip"))
    status_core <- system2(".venv/bin/pip", c("install", py_core))
    if (status_core != 0L) stop("Failed to install core Python packages into .venv.")
    try(system2(".venv/bin/pip", c("install", py_optional)), silent = TRUE)
  }
} else {
  # Ensure core packages exist inside existing .venv
  check_py <- system2(
    venv_py,
    c("-c", shQuote("import numpy, scipy, sklearn, sympy, pysindy, pandas, yaml; print('PY_OK')")),
    stdout = TRUE, stderr = TRUE
  )
  if (!any(grepl("PY_OK", check_py, fixed = TRUE))) {
    message("[INSTALL] Installing missing core Python packages into existing .venv...")
    system2(file.path(repo_root, ".venv", "bin", "pip"), c("install", py_core))
    try(system2(file.path(repo_root, ".venv", "bin", "pip"), c("install", py_optional)), silent = TRUE)
  }
}

# Generate deterministic renv.lock snapshot of active package versions
lock_pkgs <- lapply(required_r_pkgs, function(p) {
  list(Package = p, Version = as.character(packageVersion(p)), Source = "Repository")
})
names(lock_pkgs) <- required_r_pkgs
renv_lock <- list(
  R = list(Version = paste0(R.version$major, ".", R.version$minor),
           Repositories = list(list(Name = "CRAN", URL = "https://cloud.r-project.org"))),
  Packages = lock_pkgs
)
jsonlite::write_json(renv_lock, "renv.lock", pretty = TRUE, auto_unbox = TRUE)
message("[OK] Updated renv.lock.")
