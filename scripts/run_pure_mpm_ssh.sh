#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# scripts/run_pure_mpm_ssh.sh
# Audited Pure-MPM SSH Rerun Pipeline for nam-shub-01 (no SLURM required).
#
# Recommended launch on nam-shub-01 (56 cores):
#   git pull origin ICML && nohup ./scripts/run_pure_mpm_ssh.sh 44 audited_all all > run_pure_mpm.log 2>&1 &
#
# Arguments:
#   $1 NUM_WORKERS : Concurrent 1-thread R workers for real/synth (default: 44)
#   $2 SYNTH_SCOPE : Checkpoint refresh scope for Experiment 1 (synth):
#                    - audited_all (default, recommended): all 13 methods cleanly
#                      (Pure-MPM FBMS + IWSM + fixed Bag_STLSQ + deg-3 RatPoly + SR_Pool)
#                    - fbms5_sindy5: 5 FBMS methods + 5 SINDy baselines (skips IWSM to save ~8h)
#                    - fbms6: all 6 FBMS methods only
#                    - fbms5: only the 5 methods that previously had 0.025
#                    (Note: Experiment 2 'real' ALWAYS clears 'all' 2,600 tasks (~45 min)
#                     because post-split differentiation & sustained Yeast changed the data.)
#   $3 SUITE       : all (default: tuning + ablation + real + synth) |
#                    fast (tuning + ablation + real, ~1 hour total) |
#                    synth_only | real_only
# ==============================================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

NUM_WORKERS="${1:-44}"
SYNTH_SCOPE="${2:-audited_all}"
SUITE="${3:-all}"
ABL_CORES="${ABL_CORES:-6}"
export TASK_TIMEOUT="${TASK_TIMEOUT:-1800}"

export PATH="$HOME/.juliaup/bin:$HOME/.local/bin:$PATH"

# Load EasyBuild R (and Julia if available) module on nam-shub-01 (do NOT load Python module)
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

export SINDY_LEGACY_UNSTD_THRESH=0
export USE_AUTHENTIC_EMPIRICAL=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export JULIA_NUM_THREADS=1

echo "========================================================================"
echo "[PURE-MPM SSH] Starting run at $(date)"
echo "  NUM_WORKERS=$NUM_WORKERS | SYNTH_SCOPE=$SYNTH_SCOPE | SUITE=$SUITE | ABL_CORES=$ABL_CORES | TASK_TIMEOUT=$TASK_TIMEOUT"
echo "========================================================================"

echo "[STEP 0/4] Bootstrapping & verifying server R (.r_lib) + Python (.venv) + Julia dependencies..."
Rscript -e '
  lib_t <- Sys.getenv("R_LIBS_USER")
  if (!nzchar(lib_t) || !dir.exists(lib_t)) lib_t <- .libPaths()[1]
  for (pkg in c("BAS", "gauseR", "glmmTMB")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("[SETUP] Installing missing R package into ", lib_t, ": ", pkg)
      try(install.packages(pkg, lib = lib_t, repos = "https://cloud.r-project.org"), silent = TRUE)
    }
  }
'
Rscript scripts/verify_env.R

if [[ "$SUITE" == "all" || "$SUITE" == "fast" ]]; then
  echo "[STEP 1/4] Preparing real-world datasets (zero-leakage & sustained Yeast)..."
  Rscript scripts/download_real_data.R

  echo "[STEP 2/4] Re-running BGNLM calibration tuning suite (Experiment 0, Pure MPM)..."
  Rscript scripts/run_tuning.R --full --cores "$ABL_CORES"

  echo "[STEP 3/4] Re-running targeted ablation Jobs 1-4 (Pure MPM + Native PySR/PyOperon + BAS + WSINDy)..."
  ABLATION_JOBS="1,2,3,4" CORES="$ABL_CORES" ./scripts/run_targeted_ablation.sh "$ABL_CORES"
fi

# Note: In Experiments 1 & 2 (run_remote_server.sh across 44 workers), SR_Pool_Pareto & SR_Pool_MDL
# use the documented static-pool engine (PYSR_USE_JULIA=0, PYOPERON_NATIVE=0) so 44 concurrent workers
# never spawn 44 competing Julia runtimes; native Julia PySR & C++ PyOperon run in Step 3 above.
export PYSR_USE_JULIA=0
export PYOPERON_NATIVE=0

if [[ "$SUITE" == "all" || "$SUITE" == "fast" || "$SUITE" == "real_only" ]]; then
  echo "[STEP 4A] Clearing all 2,600 'real' checkpoints (post-split zero-leakage data) and running Experiment 2..."
  Rscript scripts/prepare_pure_mpm_checkpoints.R --suite real --scope all
  ./scripts/run_remote_server.sh real "$NUM_WORKERS" 1
fi

if [[ "$SUITE" == "all" || "$SUITE" == "synth_only" ]]; then
  echo "[STEP 4B] Clearing '$SYNTH_SCOPE' checkpoints for 'synth' and running Experiment 1..."
  Rscript scripts/prepare_pure_mpm_checkpoints.R --suite synth --scope "$SYNTH_SCOPE"
  ./scripts/run_remote_server.sh synth "$NUM_WORKERS" 1
fi

echo "========================================================================"
echo "[PURE-MPM SSH] ALL REQUESTED SUITES COMPLETE at $(date)."
echo "  Transfer back: results/aggregated/ results/ablation_native/ results/tables/ results/figures/ config/best_bgnlm_params.yaml"
echo "========================================================================"
