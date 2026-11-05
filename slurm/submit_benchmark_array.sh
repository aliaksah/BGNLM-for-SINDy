#!/usr/bin/env bash
#SBATCH --job-name=icml_bgnlm_synth
#SBATCH --output=logs/slurm_synth_%A_%a.out
#SBATCH --error=logs/slurm_synth_%A_%a.err
#SBATCH --array=1-1000%100
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

SUITE=${SUITE:-synth}
CHUNK_SIZE=${CHUNK_SIZE:-22}
START_ID=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
END_ID=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

TOTAL_TASKS=$(tail -n +2 "config/manifest_${SUITE}.csv" | wc -l | tr -d ' ')
if [[ "$END_ID" -gt "$TOTAL_TASKS" ]]; then
  END_ID="$TOTAL_TASKS"
fi

for TASK_ID in $(seq "$START_ID" "$END_ID"); do
    Rscript scripts/run_single_task.R --suite "$SUITE" --task_id "$TASK_ID" --cores "$SLURM_CPUS_PER_TASK"
done

