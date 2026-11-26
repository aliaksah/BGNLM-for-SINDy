#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# scripts/run_pure_mpm_ssh.sh
# End-to-end Pure-MPM SSH rerun pipeline for nam-shub-01 (no SLURM required).
# Runs fast suites first (tuning -> targeted ablation -> real -> synth),
# safely clearing old pre-Pure-MPM checkpoints once per suite so restarts resume.
#
# Usage on nam-shub-01:
#   nohup ./scripts/run_pure_mpm_ssh.sh [NUM_WORKERS=44] [SCOPE=fbms5] [SUITE=all] > run_pure_mpm.log 2>&1 &
#
# Arguments:
#   NUM_WORKERS : Concurrent xargs R workers for real/synth (default: 44 on 56-core server)
#   SCOPE       : fbms5 (default, 5 methods that had 0.025) | fbms6 (+ IWSM) |
#                 fbms5_plus_audited (5 FBMS + 7 non-FBMS audited baselines) | all
#   SUITE       : all (default: tuning + ablation + real + synth) | fast (tuning + ablation + real) |
#                 synth_only | real_only
# ==============================================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

NUM_WORKERS="${1:-44}"
SCOPE="${2:-fbms5}"
SUITE="${3:-all}"
ABL_CORES="${ABL_CORES:-8}"

# Load EasyBuild R module on nam-shub-01 if available (do NOT load Python module)
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

if [ -d "$REPO_ROOT/.r_lib" ]; then
  export R_LIBS_USER="$REPO_ROOT/.r_lib"
  export R_LIBS="$REPO_ROOT/.r_lib"
fi

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
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export JULIA_NUM_THREADS=1

echo "========================================================================"
echo "[PURE-MPM SSH] Starting run at $(date)"
echo "  NUM_WORKERS=$NUM_WORKERS | SCOPE=$SCOPE | SUITE=$SUITE | ABL_CORES=$ABL_CORES"
echo "========================================================================"

if [[ "$SUITE" == "all" || "$SUITE" == "fast" ]]; then
  echo "[STEP 1/4] Updating authentic real-world datasets..."
  Rscript scripts/download_real_data.R

  echo "[STEP 2/4] Re-running BGNLM calibration tuning suite (Experiment 0)..."
  Rscript scripts/run_tuning.R --full --cores "$ABL_CORES"

  echo "[STEP 3/4] Re-running targeted ablation Jobs 1-4 (Pure MPM + Native PySR/PyOperon)..."
  ABLATION_JOBS="1,2,3,4" CORES="$ABL_CORES" ./scripts/run_targeted_ablation.sh "$ABL_CORES"
fi

if [[ "$SUITE" == "all" || "$SUITE" == "fast" || "$SUITE" == "real_only" ]]; then
  echo "[STEP 4A] Pruning old '$SCOPE' checkpoints for 'real' and launching benchmark..."
  Rscript scripts/prepare_pure_mpm_checkpoints.R --suite real --scope "$SCOPE"
  ./scripts/run_remote_server.sh real "$NUM_WORKERS" 1
fi

if [[ "$SUITE" == "all" || "$SUITE" == "synth_only" ]]; then
  echo "[STEP 4B] Pruning old '$SCOPE' checkpoints for 'synth' and launching benchmark..."
  Rscript scripts/prepare_pure_mpm_checkpoints.R --suite synth --scope "$SCOPE"
  ./scripts/run_remote_server.sh synth "$NUM_WORKERS" 1
fi

echo "========================================================================"
echo "[PURE-MPM SSH] Finished at $(date). Aggregated results in results/aggregated/ and results/ablation_native/."
echo "========================================================================"
