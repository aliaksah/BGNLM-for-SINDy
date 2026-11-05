#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/run_remote_server.sh [SUITE: synth|real|both] [NUM_WORKERS] [THREADS_PER_WORKER] [START_TASK] [END_TASK]
# Example (run via nohup over SSH):
#   nohup ./scripts/run_remote_server.sh both 16 2 > run_benchmarks.log 2>&1 &
SUITE=${1:-both}
NUM_WORKERS=${2:-8}
THREADS_PER_WORKER=${3:-4}
START_TASK_USER=${4:-}
END_TASK_USER=${5:-}

# Load environment modules if available (e.g. on HPC / shared servers)
if command -v module >/dev/null 2>&1; then
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
  module load Python 2>/dev/null || true
elif [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
  module load Python 2>/dev/null || true
fi

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

mkdir -p results/checkpoints_synth results/checkpoints_real logs

echo "[INFO] Ensuring real-world datasets are present..."
Rscript scripts/download_real_data.R

if [[ ! -f "config/manifest_synth.csv" || ! -f "config/manifest_real.csv" ]]; then
  echo "[INFO] Generating default manifests..."
  Rscript scripts/generate_manifest.R
fi

run_suite() {
  local s=$1
  local manifest="config/manifest_${s}.csv"
  local total_tasks
  total_tasks=$(tail -n +2 "$manifest" | wc -l | tr -d ' ')

  local start_task=${START_TASK_USER:-1}
  local end_task=${END_TASK_USER:-$total_tasks}

  if (( end_task > total_tasks )); then
    end_task=$total_tasks
  fi

  echo "[INFO] Launching tasks $start_task to $end_task for suite '$s' across $NUM_WORKERS parallel workers ($THREADS_PER_WORKER threads each)..."
  seq "$start_task" "$end_task" | xargs -n 1 -P "$NUM_WORKERS" -I {} bash -c "
    Rscript scripts/run_single_task.R --suite $s --task_id {} --cores $THREADS_PER_WORKER > logs/${s}_task_{}.log 2>&1
  "
}

if [[ "$SUITE" == "synth" || "$SUITE" == "both" ]]; then
  run_suite "synth"
fi

if [[ "$SUITE" == "real" || "$SUITE" == "both" ]]; then
  run_suite "real"
fi

echo "[INFO] Aggregating checkpoints and running mixed-effects evaluation..."
Rscript scripts/aggregate_results.R
Rscript scripts/run_mixed_models.R
echo "[INFO] Pipeline complete."

