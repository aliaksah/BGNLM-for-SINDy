#!/usr/bin/env bash
set -euo pipefail

# Targeted server ablation: native PySR/PyOperon + standardized UQ_SINDy + state_noise + authentic 1900-1920 Lynx-Hare
CORES="${1:-6}"

export PYSR_USE_JULIA=1
export PYOPERON_NATIVE=1
export SINDY_LEGACY_UNSTD_THRESH=0
export USE_AUTHENTIC_EMPIRICAL=1
export CORES="${CORES}"

if [ -f ".venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source .venv/bin/activate
fi

echo "[ABLATION] Starting targeted server ablation with CORES=${CORES}..."
Rscript scripts/run_targeted_ablation.R
echo "[ABLATION] Complete. Outputs in results/ablation_native/."
