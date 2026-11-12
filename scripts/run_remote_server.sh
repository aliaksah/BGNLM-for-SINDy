#!/usr/bin/env bash
set -euo pipefail

# Usage: ./scripts/run_remote_server.sh [SUITE: synth|real|both] [NUM_WORKERS] [THREADS_PER_WORKER] [START_TASK] [END_TASK]
# Recommended on nam-shub-01 (56 cores) via nohup:
#   nohup ./scripts/run_remote_server.sh both 40 1 > run_benchmarks.log 2>&1 &

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

SUITE=${1:-both}
NUM_WORKERS=${2:-36}
THREADS_PER_WORKER=${3:-1}
START_TASK_USER=${4:-}
END_TASK_USER=${5:-}

# Load R module ONLY if available (do NOT load Python module as it can conflict with .venv)
if command -v module >/dev/null 2>&1; then
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
elif [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /etc/profile.d/modules.sh
  module load R/4.5.2-gfbf-2025b 2>/dev/null || true
fi

if ! command -v Rscript >/dev/null 2>&1; then
  echo "[ERROR] Rscript not found! Run ./scripts/setup_server.sh or load R module first."
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

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export JULIA_NUM_THREADS=1

mkdir -p results/checkpoints_synth results/checkpoints_real logs

echo "[INFO] Running pre-flight environment verification..."
Rscript scripts/verify_env.R

echo "[INFO] Ensuring real-world datasets are present..."
Rscript scripts/download_real_data.R

if [[ ! -f "config/manifest_synth.csv" || ! -f "config/manifest_real.csv" ]]; then
  echo "[INFO] Generating default manifests..."
  Rscript scripts/generate_manifest.R
fi

# Resolve OS timeout command if available (GNU coreutils on Linux)
TIMEOUT_CMD=""
if command -v timeout >/dev/null 2>&1; then
  TIMEOUT_CMD="timeout -k 15 360"
fi
export TIMEOUT_CMD

run_suite() {
  local s=$1
  local manifest="config/manifest_${s}.csv"
  local ckpt_dir="results/checkpoints_${s}"
  local total_tasks
  total_tasks=$(tail -n +2 "$manifest" | wc -l | tr -d ' ')

  local start_task=${START_TASK_USER:-1}
  local end_task=${END_TASK_USER:-$total_tasks}
  if (( end_task > total_tasks )); then
    end_task=$total_tasks
  fi

  # Emit ONLY missing task_ids in a single fast pass (< 0.5s) to avoid spawning R for completed tasks
  local pending_file="logs/.pending_${s}_tasks.txt"
  Rscript -e "
    dir_p <- '$ckpt_dir'
    existing <- list.files(dir_p, pattern = '^task_[0-9]+[.]rds$')
    done_ids <- if (length(existing) > 0L) as.integer(gsub('[^0-9]', '', existing)) else integer(0)
    all_ids <- seq.int(${start_task}L, ${end_task}L)
    todo <- setdiff(all_ids, done_ids)
    writeLines(as.character(todo), '$pending_file')
    cat(sprintf('[INFO] Suite %s: %d total in range [%d..%d], %d already done, %d pending.',
                '$s', length(all_ids), ${start_task}L, ${end_task}L, length(all_ids) - length(todo), length(todo)), '\n')
  "

  local n_pending
  n_pending=$(wc -l < "$pending_file" | tr -d ' ')
  if (( n_pending == 0 )); then
    echo "[INFO] Suite '$s' is already 100% complete!"
    return 0
  fi

  local fail_log="logs/${s}_failures.log"
  : > "$fail_log"

  echo "[INFO] Launching $n_pending pending tasks for suite '$s' across $NUM_WORKERS workers ($THREADS_PER_WORKER thread(s) each)..."
  export SUITE_CUR="$s"
  export THREADS_PER_WORKER
  export FAIL_LOG="$fail_log"

  run_one_task() {
    local tid="$1"
    local tlog="logs/${SUITE_CUR}_task_${tid}.log"
    if ${TIMEOUT_CMD} Rscript scripts/run_single_task.R --suite "${SUITE_CUR}" --task_id "${tid}" --cores "${THREADS_PER_WORKER}" --timeout 300 > "${tlog}" 2>&1; then
      rm -f "${tlog}"
    else
      local rc=$?
      echo "[FAIL] suite=${SUITE_CUR} task_id=${tid} exit=${rc} (see ${tlog})" >> "${FAIL_LOG}"
    fi
  }
  export -f run_one_task

  xargs -n 1 -P "$NUM_WORKERS" bash -c 'run_one_task "$1"' _ < "$pending_file"

  rm -f "$pending_file"
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
