#!/usr/bin/env bash
# Data requirements: pull every summary statistic the analysis needs from OpenGWAS into data/raw/.
#
#   export OPENGWAS_JWT="eyJ..."                  # free token: https://api.opengwas.io/profile/
#   export KEATON_XLSX=/path/to/keaton_supp.xlsx  # optional: Keaton et al. 2024 supplementary tables
#   bash fetch_data.sh
#
# Needs internet access (on Kaggle: Settings -> Internet -> On). Safe to re-run: each study's extract
# (and every API call) is saved in data/raw/cache/, so an interrupted run resumes where it stopped.
# Options: FETCH_WORKERS=3 (studies at once), KEEP_FILES_DIR=/path (also keep the full .vcf.gz files),
#          FETCH_SOURCE=api (old per-variant look-ups).
set -euo pipefail
cd "$(dirname "$0")"

if [[ -z "${OPENGWAS_JWT:-}" ]]; then
  echo "ERROR: OPENGWAS_JWT is not set. Get a token at https://api.opengwas.io/profile/" >&2
  exit 1
fi
echo "Checking OpenGWAS is reachable..."
curl -sf -m 30 -H "Authorization: Bearer ${OPENGWAS_JWT}" https://api.opengwas.io/api/status > /dev/null \
  || { echo "ERROR: cannot reach api.opengwas.io (internet off, or token invalid)" >&2; exit 1; }

ARGS=()
if [[ -n "${KEATON_XLSX:-}" ]]; then
  if [[ -f "${KEATON_XLSX}" ]]; then ARGS+=(--keaton "${KEATON_XLSX}"); else echo "WARNING: KEATON_XLSX=${KEATON_XLSX} not found; continuing without it"; fi
fi

# FETCH_SOURCE=files (default: stream whole study files, ~1-2 h) or api (per-variant look-ups, LD proxies, ~2 min per 64)
ARGS+=(--source "${FETCH_SOURCE:-files}" --workers "${FETCH_WORKERS:-3}")
if [[ -n "${KEEP_FILES_DIR:-}" ]]; then ARGS+=(--keep-files "${KEEP_FILES_DIR}"); fi
python scripts/00_fetch_opengwas.py "${ARGS[@]+"${ARGS[@]}"}"

for f in manifest.csv associations.csv.gz instruments.csv; do
  [[ -s "data/raw/$f" ]] || { echo "ERROR: data/raw/$f missing or empty" >&2; exit 1; }
done
echo
echo "Studies resolved (check these before trusting the results):"
python - <<'PY'
import pandas as pd
m = pd.read_csv("data/raw/manifest.csv")
cols = [c for c in ["key", "role", "id", "trait", "author", "year", "sample_size", "ncase", "unit", "how"] if c in m]
print(m[cols].to_string(index=False))
missing = m[m["id"].isna()]["key"].tolist()
if missing:
    print("\nNOT FOUND (analysis runs without them):", missing)
import os
if os.path.exists("data/raw/file_extract_report.csv"):
    r = pd.read_csv("data/raw/file_extract_report.csv")
    print("\nVariants missing per study (no LD proxies in file mode), worst first:")
    print(r.head(8).to_string(index=False))
if os.path.exists("data/raw/rejected_studies.csv"):
    print("\nStudies that could not be fetched:")
    print(pd.read_csv("data/raw/rejected_studies.csv").to_string(index=False))
PY
python - <<'PY'
import zipfile, pathlib
root = pathlib.Path("data")
with zipfile.ZipFile(root / "raw_fetched.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for f in (root / "raw").rglob("*"):
        if f.is_file() and "cache" not in f.parts:
            z.write(f, f.relative_to(root))
PY
echo
echo "Done. data/raw/ is ready; a copy without the cache is in data/raw_fetched.zip"
