# ==============================================================================
# R/00_env_setup.R
# Enforces Project-Local R Library (.r_lib) & Python Virtualenv (.venv)
# Guarantees deterministic, self-contained execution on shared Linux/HPC servers
# ==============================================================================

setup_project_env <- function() {
  # Locate repository root (works from repo root or tests/testthat)
  candidates <- c(".", "..", "../..")
  repo_root <- "."
  for (cand in candidates) {
    if (file.exists(file.path(cand, "DESCRIPTION")) && dir.exists(file.path(cand, "config"))) {
      repo_root <- normalizePath(cand, mustWork = TRUE)
      break
    }
  }

  local_r_lib <- file.path(repo_root, ".r_lib")
  if (dir.exists(local_r_lib)) {
    local_r_lib <- normalizePath(local_r_lib, mustWork = TRUE)
    # Place project-local .r_lib strictly first, followed by base .Library and existing paths
    .libPaths(unique(c(local_r_lib, .Library, .libPaths())))
    Sys.setenv(R_LIBS_USER = local_r_lib, R_LIBS = local_r_lib)
  }

  venv_py <- file.path(repo_root, ".venv", "bin", "python")
  if (file.exists(venv_py)) {
    venv_py_dir <- normalizePath(dirname(venv_py), mustWork = TRUE)
    venv_py_abs <- file.path(venv_py_dir, "python")
    Sys.setenv(RETICULATE_PYTHON = venv_py_abs)
    cur_path <- Sys.getenv("PATH")
    if (!startsWith(cur_path, venv_py_dir)) {
      Sys.setenv(PATH = paste(venv_py_dir, cur_path, sep = .Platform$path.sep))
    }
  }

  invisible(list(repo_root = repo_root, local_r_lib = local_r_lib))
}

setup_project_env()

