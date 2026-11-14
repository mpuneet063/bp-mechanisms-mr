#!/usr/bin/env bash
# Install Python dependencies (Kaggle already has most of them; this is quick).
set -euo pipefail
cd "$(dirname "$0")"
python -m pip install -q -r requirements.txt
python -c "import numpy, pandas, scipy, sklearn, matplotlib, requests, openpyxl; print('dependencies OK')"
