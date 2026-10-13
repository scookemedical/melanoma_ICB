#!/usr/bin/env bash
# Acquire the two external inputs. Both come from GitHub because this session's
# egress policy blocks ftp.ncbi.nlm.nih.gov, ebi.ac.uk, zenodo.org and
# cbioportal.org (see docs/METHODS.md).
set -euo pipefail
DEST="${1:?usage: 01_fetch_data.sh <download-dir>}"
mkdir -p "$DEST/humangem" "$DEST/riaz"

# 1. Human-GEM v2.1.1 -- genome-scale human metabolic reconstruction
BASE=https://raw.githubusercontent.com/SysBioChalmers/Human-GEM/main
for f in model/Human-GEM.xml model/genes.tsv model/reactions.tsv model/metabolites.tsv; do
  curl -fsSL --retry 4 --retry-delay 2 -o "$DEST/humangem/$(basename "$f")" "$BASE/$f"
done

# 2. Riaz et al., Cell 2017 (GSE91061) -- nivolumab-treated metastatic melanoma,
#    paired pre/on-treatment bulk RNA-seq, from the authors' own analysis repo.
if [ ! -d "$DEST/riaz/.git" ]; then
  git clone --depth 1 https://github.com/riazn/bms038_analysis.git "$DEST/riaz"
fi
echo "inputs ready in $DEST"
