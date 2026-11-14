#!/usr/bin/env bash
# Sanity check before touching real data: unit tests, then the full pipeline on simulated data with known truth.
#   bash run_simulation_check.sh          # full bootstrap counts (~10-20 min on Kaggle CPU)
#   bash run_simulation_check.sh --fast   # ~1 min
set -euo pipefail
cd "$(dirname "$0")"
python -m pytest -q tests
python -m bpmr.simulate data/sim_raw
OUT="results_sim"; [[ "${1:-}" == "--fast" ]] && OUT="results_sim_fast"
bash run_experiments.sh --raw data/sim_raw --out "$OUT" "$@"
echo; cat "$OUT/sim_validation.md"
