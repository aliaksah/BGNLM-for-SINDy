#!/usr/bin/env bash
set -euo pipefail

# Targeted server ablation: native PySR/PyOperon + standardized UQ_SINDy + state_noise + authentic 1900-1920 Lynx-Hare
# Recommended via nohup:
#   nohup ./scripts/run_targeted_ablation.sh 6 > run_targeted_ablation.log 2>&1 &

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

CORES="${1:-6}"

# Load EasyBuild R module on nam-shub-01 if available
if command -v module >/dev/null 2>&1; then
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
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

# Activate project-local Python .venv if present
if [ -f "$REPO_ROOT/.venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "$REPO_ROOT/.venv/bin/activate"
  export VIRTUAL_ENV="$REPO_ROOT/.venv"
  export PATH="$REPO_ROOT/.venv/bin:$PATH"
fi

export PYSR_USE_JULIA=1
export PYOPERON_NATIVE=1
export SINDY_LEGACY_UNSTD_THRESH=0
export USE_AUTHENTIC_EMPIRICAL=1
export CORES="${CORES}"

echo "[ABLATION] Starting targeted server ablation with CORES=${CORES}..."
Rscript scripts/run_targeted_ablation.R
echo "[ABLATION] Complete. Outputs in results/ablation_native/."
