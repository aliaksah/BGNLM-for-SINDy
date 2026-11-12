#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# scripts/setup_server.sh
# One-shot bootstrap for remote Linux servers (e.g., nam-shub-01):
# Creates self-contained project-local R (.r_lib/) and Python (.venv/) environments
# ==============================================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "[SETUP] Repository root: $REPO_ROOT"

# Load R module on EasyBuild/Lmod clusters if Rscript isn't already available or module exists
if command -v module >/dev/null 2>&1; then
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found! Please run 'module load R/4.5.2-gfbf-2025b' (or equivalent) first."
  exit 1
fi

echo "[SETUP] Using Rscript: $(command -v Rscript)"
mkdir -p "$REPO_ROOT/.r_lib" results/checkpoints_synth results/checkpoints_real logs

export R_LIBS_USER="$REPO_ROOT/.r_lib"
export R_LIBS="$REPO_ROOT/.r_lib"
export R_MAKEVARS_USER=/dev/null

# If .venv exists, activate it first so install_dependencies.R uses it
if [ -f "$REPO_ROOT/.venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.venv/bin/activate"
fi

echo "[SETUP] Installing R packages into $REPO_ROOT/.r_lib and verifying .venv..."
Rscript scripts/install_dependencies.R

echo "[SETUP] Running pre-flight environment verification..."
Rscript scripts/verify_env.R

echo "[SETUP] Preparing real-world datasets..."
Rscript scripts/download_real_data.R

echo "[SETUP] Generating manifests if missing..."
if [[ ! -f "config/manifest_synth.csv" || ! -f "config/manifest_real.csv" ]]; then
  Rscript scripts/generate_manifest.R
fi

echo "[SETUP] Running unit test suite..."
Rscript -e 'testthat::test_dir("tests/testthat")'

echo "[OK] Server setup complete! You can now launch benchmarks with:"
echo "     nohup ./scripts/run_remote_server.sh both 40 1 > run_benchmarks.log 2>&1 &"

