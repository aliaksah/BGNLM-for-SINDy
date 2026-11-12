#!/usr/bin/env Rscript
# ==============================================================================
# scripts/verify_env.R
# Fast Pre-Flight Verification of Project-Local R (.r_lib) & Python (.venv)
# ==============================================================================

if (file.exists("R/00_env_setup.R")) source("R/00_env_setup.R")

message("[VERIFY] Active R version: ", R.version.string)
message("[VERIFY] Active .libPaths():\n  - ", paste(.libPaths(), collapse = "\n  - "))

required_r_pkgs <- c(
  "FBMS", "deSolve", "glmnet", "signal", "data.table", "yaml", "jsonlite",
  "lme4", "lmerTest", "glmmTMB", "emmeans", "broom.mixed",
  "ggplot2", "dplyr", "tidyr", "reticulate", "R.utils", "testthat"
)

for (pkg in required_r_pkgs) {
  ok <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
    TRUE
  }, error = function(e) {
    message("[FAIL] Cannot load R package '", pkg, "': ", e$message)
    FALSE
  })
  if (!ok) stop("Pre-flight R check failed on package: ", pkg)
}

# Verify Python .venv bridge
source("R/04b_sr_baselines.R")
py_bin <- find_python_bin()
message("[VERIFY] Python interpreter: ", py_bin)

py_chk <- system2(
  py_bin,
  c("-c", shQuote("import numpy, scipy, sklearn, sympy, pysindy, pandas, yaml; print('PYTHON_ENV_OK')")),
  stdout = TRUE, stderr = TRUE
)
if (!any(grepl("PYTHON_ENV_OK", py_chk, fixed = TRUE))) {
  stop("Pre-flight Python check failed! Output:\n", paste(py_chk, collapse = "\n"))
}

# Quick sanity check on FBMS + symbolic parser
source("R/03_fbms_wrappers.R")
t_dicts <- register_fbms_transforms()
deep_test <- extract_canonical_atoms("p3(1+1*x1+1*x2+1*p3(1+1*x1+1*x2+1*p2(x1)))")
stopifnot(length(deep_test) >= 1L)

message("[OK] Pre-flight verification PASSED (R + FBMS ", as.character(packageVersion("FBMS")), " + Python .venv ready).")

