#!/usr/bin/env Rscript
# ==============================================================================
# scripts/install_dependencies.R
# Installs Required R Packages & Configures Python Virtual Environment (.venv)
# ==============================================================================

# Ensure a writable user library directory exists (standard on shared Linux/HPC clusters)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!nzchar(user_lib)) {
  user_lib <- path.expand(file.path("~", "R", paste0(R.version$platform, "-library"), paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1])))
}
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(c(user_lib, .libPaths()))

required_r_pkgs <- c(
  "FBMS", "deSolve", "glmnet", "signal", "data.table", "yaml", "jsonlite",
  "lme4", "lmerTest", "glmmTMB", "emmeans", "broom.mixed",
  "ggplot2", "dplyr", "tidyr", "reticulate", "R.utils", "testthat"
)

for (pkg in required_r_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("[INSTALL] Installing R package: ", pkg, " into ", user_lib)
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org")
  }
}

message("[OK] All required R packages are installed (FBMS version: ",
        as.character(packageVersion("FBMS")), ").")

# Set up Python environment if not already present
if (!file.exists(".venv/bin/python")) {
  uv_bin <- Sys.which("uv")
  if (!nzchar(uv_bin) && file.exists( path.expand("~/.local/bin/uv") )) {
    uv_bin <- path.expand("~/.local/bin/uv")
  }
  if (nzchar(uv_bin)) {
    message("[INSTALL] Creating .venv via uv...")
    system2(uv_bin, c("venv", ".venv", "--python", "3.11"))
    system2(uv_bin, c("pip", "install", "--python", ".venv/bin/python",
                      "numpy", "scipy", "scikit-learn", "sympy",
                      "pysindy", "pysr", "pyoperon", "pandas", "pyyaml"))
  } else {
    py3 <- Sys.which("python3")
    if (nzchar(py3)) {
      system2(py3, c("-m", "venv", ".venv"))
      system2(".venv/bin/pip", c("install", "numpy", "scipy", "scikit-learn", "sympy",
                                 "pysindy", "pysr", "pyoperon", "pandas", "pyyaml"))
    }
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

