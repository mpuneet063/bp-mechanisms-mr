#!/usr/bin/env bash
# One-shot Kaggle entry point: setup -> tests + simulation check -> fetch real data -> all experiments -> zip.
#
# In a Kaggle notebook (Internet ON, secret OPENGWAS_JWT added under Add-ons -> Secrets):
#   import os; from kaggle_secrets import UserSecretsClient
#   os.environ["OPENGWAS_JWT"] = UserSecretsClient().get_secret("OPENGWAS_JWT")
#   os.environ["KEATON_XLSX"] = "/kaggle/input/<your-dataset>/keaton_supp.xlsx"   # optional
#   !bash /kaggle/input/<your-dataset>/bp-mechanisms-mr/kaggle_run.sh
#
# Set SKIP_SIM=1 to skip the simulation check, FAST=1 for small bootstrap counts.
# Safe to re-run after a crash: the fetch resumes from data/raw/cache and finished analysis steps are skipped.
set -euo pipefail
SRC="$(cd "$(dirname "$0")" && pwd)"
# Kaggle: /kaggle/input is read-only, so work on a copy in /kaggle/working.
# Anywhere else (Colab with the repo on Google Drive, a laptop): run in place, so progress persists.
if [[ "$SRC" == /kaggle/input/* ]]; then DEFAULT_WORK=/kaggle/working/bp-mechanisms-mr; else DEFAULT_WORK="$SRC"; fi
WORK="${WORK:-$DEFAULT_WORK}"
if [[ "$SRC" != "$WORK" ]]; then
  mkdir -p "$WORK"; cp -rn "$SRC/." "$WORK/"   # -n: never overwrite progress already in the working copy
fi
cd "$WORK"
chmod +x ./*.sh
bash setup.sh
FASTFLAG=(); [[ "${FAST:-0}" == "1" ]] && FASTFLAG=(--fast)
if [[ "${SKIP_SIM:-0}" != "1" && ! -f results_sim_fast/sim_validation.md ]]; then
  bash run_simulation_check.sh --fast
fi
bash fetch_data.sh
bash run_experiments.sh "${FASTFLAG[@]+"${FASTFLAG[@]}"}"
python - <<'PY'
import zipfile, pathlib, os
work = pathlib.Path(os.getcwd())
dest = work.parent / "results_bundle.zip"
with zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as z:
    for f in list((work / "results").rglob("*")) + [work / "data/raw" / n for n in
                  ("manifest.csv", "manifest_alternatives.csv", "fetch_log.txt", "instruments.csv")]:
        if f.is_file() and f.suffix != ".npy":
            z.write(f, f.relative_to(work.parent))
print("wrote", dest)
PY
echo "All done. Download /kaggle/working/results_bundle.zip (results + study manifest)."
