#!/usr/bin/env bash
#SBATCH --job-name=icml_bgnlm_tune
#SBATCH --output=logs/slurm_tune_%j.out
#SBATCH --error=logs/slurm_tune_%j.err
#SBATCH --cpus-per-task=10
#SBATCH --mem=16G
#SBATCH --time=04:00:00

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

mkdir -p logs results/aggregated results/tables
Rscript scripts/run_tuning.R --full --cores "${SLURM_CPUS_PER_TASK:-10}"

