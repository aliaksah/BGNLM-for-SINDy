#!/usr/bin/env bash
set -euo pipefail

# Targeted server ablation: Pure-MPM BGNLM + native PySR/PyOperon + standardized Bag_STLSQ
# + BAS Bayesian SINDy + Galerkin WSINDy + authentic zero-leakage empirical datasets
# Recommended via nohup:
#   nohup ./scripts/run_targeted_ablation.sh 6 > run_targeted_ablation.log 2>&1 &

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CORES="${1:-6}"

export PATH="$HOME/.juliaup/bin:$HOME/.local/bin:$PATH"

# Load EasyBuild R (and Julia if available) module on nam-shub-01
if command -v module >/dev/null 2>&1; then
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
  module load Julia 2>/dev/null || true
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
  module load Julia 2>/dev/null || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found! Load the R module or run ./scripts/setup_server.sh first."
  exit 1
fi

# Enforce project-local R library (.r_lib) if present
if [ -d "$REPO_ROOT/.r_lib" ]; then
  export R_LIBS_USER="$REPO_ROOT/.r_lib"
  export R_LIBS="$REPO_ROOT/.r_lib"
fi

# Ensure BAS (and optional gauseR) exist in R library for Job 2 & Job 4
Rscript -e '
  lib_t <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib_t) || !dir.exists(lib_t)) lib_t <- .libPaths()[1]
  for (pkg in c("BAS", "gauseR")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[ABLATION SETUP] Installing missing R package: ", pkg)
      try(install.packages(pkg, lib = lib_t, repos = "https://cloud.r-project.org"), silent = TRUE)
    }
  }
  stopifnot(requireNamespace("BAS", quietly = TRUE))
'

# Activate project-local Python .venv if present
if [ -f "$REPO_ROOT/.venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.venv/bin/activate"
  export VIRTUAL_ENV="$REPO_ROOT/.venv"
  export PATH="$REPO_ROOT/.venv/bin:$PATH"
fi

# Ensure pyoperon and pysr are installed inside .venv
PY_BIN="${VIRTUAL_ENV:+$VIRTUAL_ENV/bin/python}"
PY_BIN="${PY_BIN:-python3}"
if ! "$PY_BIN" -c "import pyoperon.sklearn, pysr" >/dev/null 2>&1; then
  echo "[ABLATION SETUP] Installing pyoperon & pysr into .venv..."
  if command -v uv >/dev/null 2>&1 && [ -n "${VIRTUAL_ENV:-}" ]; then
    uv pip install --python "$PY_BIN" pyoperon pysr || true
  elif [ -n "${VIRTUAL_ENV:-}" ]; then
    "$VIRTUAL_ENV/bin/pip" install pyoperon pysr || true
  fi
fi

# Ensure Julia is available for native PySR
if ! command -v julia >/dev/null 2>&1; then
  echo "[ABLATION SETUP] 'julia' not found on PATH; bootstrapping standalone juliaup in ~/.juliaup..."
  curl -fsSL https://install.julialang.org | sh -s -- -y --add-to-path=no || true
  export PATH="$HOME/.juliaup/bin:$PATH"
fi

if command -v julia >/dev/null 2>&1; then
  "$PY_BIN" -c "import pysr; pysr.install(quiet=True)" >/dev/null 2>&1 || true
fi

export PYSR_USE_JULIA=1
export PYOPERON_NATIVE=1
export SINDY_LEGACY_UNSTD_THRESH=0
export USE_AUTHENTIC_EMPIRICAL=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
export CORES="${CORES}"

echo "[ABLATION PREFLIGHT] Checking native PyOperon and PySR/Julia bridges..."
"$PY_BIN" -c '
import shutil
try:
    from pyoperon.sklearn import SymbolicRegressor
    print("  [OK] pyoperon native C++ backend ready.")
except Exception as e:
    print(f"  [WARN] pyoperon unavailable: {e}")
try:
    import pysr
    print(f"  [OK] pysr package ready (julia={shutil.which(\"julia\")}).")
except Exception as e:
    print(f"  [WARN] pysr unavailable: {e}")
'

echo "[ABLATION] Starting targeted server ablation with CORES=${CORES}..."
Rscript scripts/run_targeted_ablation.R
echo "[ABLATION] Complete. Outputs in results/ablation_native/."
