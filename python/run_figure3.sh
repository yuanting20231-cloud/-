#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python - <<'PY'
from pathlib import Path
import base64
src = Path('data/figure3_family_counts.tsv.gz.b64')
dst = Path('data/figure3_family_counts.tsv.gz')
dst.write_bytes(base64.b64decode(src.read_text(encoding='utf-8').strip()))
PY
python - <<'PY'
import pandas as pd
pd.read_csv('data/figure3_family_counts.tsv.gz', sep='\t', compression='gzip').to_csv(
    'data/figure3_family_counts.tsv', sep='\t', index=False
)
PY
python python/run_figure3.py
