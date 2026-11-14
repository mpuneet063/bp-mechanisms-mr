#!/usr/bin/env bash
# Data requirements: pull every summary statistic the analysis needs from OpenGWAS into data/raw/.
#
#   export OPENGWAS_JWT="eyJ..."                  # free token: https://api.opengwas.io/profile/
#   export KEATON_XLSX=/path/to/keaton_supp.xlsx  # optional: Keaton et al. 2024 supplementary tables
#   bash fetch_data.sh
#
# Needs internet access (on Kaggle: Settings -> Internet -> On). Safe to re-run: API calls are cached
# in data/raw/cache/, so an interrupted run resumes where it stopped.
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
