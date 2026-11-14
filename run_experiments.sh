#!/usr/bin/env bash
# Run every experiment in order. Usage:
#   bash run_experiments.sh                       # real data in data/raw -> results/
#   bash run_experiments.sh --fast                # same, small bootstrap counts (quick smoke test)
#   bash run_experiments.sh --raw data/sim_raw --out results_sim
# Extra flags (--primary bnmf, --k 6, --profile raw, --phi 0.25, --seed N) are passed to every step.
set -euo pipefail
cd "$(dirname "$0")"

RAW="data/raw"; OUT="results"; EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --raw) RAW="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done
mkdir -p "$OUT"
LOG="$OUT/run.log"
echo "run started $(date '+%F %T')  raw=$RAW out=$OUT extra=${EXTRA[*]-}" | tee "$LOG"

for step in 01_prepare 02_cluster 03_stability 04_mr 05_compare 06_mvmr 07_report; do
  echo "=== $step  $(date '+%T')" | tee -a "$LOG"
  python "scripts/$step.py" --raw "$RAW" --out "$OUT" "${EXTRA[@]+"${EXTRA[@]}"}" 2>&1 | tee -a "$LOG"
done
if [[ -f "$RAW/truth.csv" ]]; then
  python scripts/sim_validate.py --raw "$RAW" --out "$OUT" "${EXTRA[@]+"${EXTRA[@]}"}" 2>&1 | tee -a "$LOG"
fi
echo "run finished $(date '+%F %T')" | tee -a "$LOG"
echo "Summary: $OUT/07_report/summary.md   Figures: $OUT/07_report/"
