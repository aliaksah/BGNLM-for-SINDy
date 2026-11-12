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
n_cpus <- min(8L, max(1L, n_cores - 1L))
options(Ncpus = n_cpus)
message("[ENV] Using Ncpus = ", n_cpus, " for parallel source package compilation.")

# CRAN packages needed (remotes added to install from GitHub branches)
required_cran_pkgs <- c(
  "remotes", "deSolve", "glmnet", "signal", "data.table", "yaml", "jsonlite",
  "lme4", "lmerTest", "glmmTMB", "emmeans", "broom.mixed",
  "ggplot2", "dplyr", "tidyr", "reticulate", "R.utils", "testthat"
)

is_pkg_usable <- function(pkg, lib_loc) {
  # Check that package exists in target_lib (or base .Library) AND loads cleanly
  found <- nzchar(system.file(package = pkg, lib.loc = c(lib_loc, .Library)))
  if (!found) return(FALSE)
  tryCatch({
    suppressPackageStartupMessages(loadNamespace(pkg, lib.loc = c(lib_loc, .Library)))
    TRUE
  }, error = function(e) FALSE)
}

missing_cran <- required_cran_pkgs[!vapply(required_cran_pkgs, is_pkg_usable, logical(1), lib_loc = target_lib)]

if (length(missing_cran) > 0L) {
  message("[INSTALL] Installing missing CRAN packages (", length(missing_cran), "): ",
          paste(missing_cran, collapse = ", "))
  install.packages(
    missing_cran,
    lib = target_lib,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE,
    Ncpus = n_cpus
  )
}

# Install / update FBMS specifically from GitHub branch: jonlachmann/FBMS@JCS
message("[INSTALL] Installing FBMS from GitHub: jonlachmann/FBMS@JCS into ", target_lib, "...")
remotes::install_github(
  "jonlachmann/FBMS@JCS",
  lib = target_lib,
  upgrade = "never",
  force = TRUE,
  quiet = FALSE
)

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
