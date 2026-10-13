#!/usr/bin/env bash
# End-to-end run. WORK defaults to a scratch dir because the flux sweep and the
# external inputs are large and do not belong in the repository.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$HERE")"
WORK="${WORK:-${TMPDIR:-/tmp}/gln_rank}"
DL="$WORK/dl"
OUT="$REPO/results"
TIMEPOINT="${TIMEPOINT:-Pre}"
DOSES="${DOSES:-1.0,1.25,1.5,2.0,2.5}"
VMAX="${VMAX:-0.3}"      # calibrated; see docs/METHODS.md
FLOOR="${FLOOR:-0.005}"
PROCS="${PROCS:-4}"
mkdir -p "$WORK" "$OUT"

bash "$HERE/01_fetch_data.sh" "$DL"

python3 "$HERE/00_load_model.py" "$DL/humangem/Human-GEM.xml" "$WORK/humangem.pkl"

python3 "$HERE/02_prepare_expression.py" \
  --counts   "$DL/riaz/data/CountData.BMS038.txt" \
  --clinical "$DL/riaz/data/bms038_clinical_data.csv" \
  --gem-genes "$DL/humangem/genes.tsv" \
  --out-expr "$WORK/expr_logcpm_ensg.tsv" \
  --out-samples "$WORK/samples.tsv"

python3 "$HERE/03_eflux_sweep.py" \
  --model "$WORK/humangem.pkl" \
  --expr  "$WORK/expr_logcpm_ensg.tsv" \
  --samples "$WORK/samples.tsv" \
  --panel "$REPO/config/metabolite_panel.tsv" \
  --out "$WORK/sweep_${TIMEPOINT}.tsv" \
  --timepoint "$TIMEPOINT" --doses "$DOSES" --vmax "$VMAX" \
  --floor-frac "$FLOOR" --procs "$PROCS"

python3 "$HERE/04_rank_metabolites.py" \
  --sweep "$WORK/sweep_${TIMEPOINT}.tsv" \
  --panel "$REPO/config/metabolite_panel.tsv" \
  --out-ranking  "$OUT/metabolite_ranking_${TIMEPOINT}.tsv" \
  --out-stability "$OUT/rank_stability_${TIMEPOINT}.tsv"

python3 "$HERE/05_validate_icb.py" \
  --sweep "$WORK/sweep_${TIMEPOINT}.tsv" \
  --samples "$WORK/samples.tsv" \
  --panel "$REPO/config/metabolite_panel.tsv" \
  --ranking "$OUT/metabolite_ranking_${TIMEPOINT}.tsv" \
  --out "$OUT/icb_association_${TIMEPOINT}.tsv"

echo "done -- results in $OUT"
