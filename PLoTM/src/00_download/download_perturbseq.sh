#!/bin/bash
# download_perturbseq.sh — fetch the CRISPRi Perturb-seq dataset that anchors the pipeline.
#
# Default target = the genome-scale primary human CD4+ T-cell Perturb-seq (CZI virtualcellmodels),
# the dataset used for the RA reference run. For a DIFFERENT disease/cell type, point PERTURB_URL /
# PERTURB_OUT at that dataset's pseudobulk h5ad (see config/disease.example.yaml -> perturbseq.url).
#
# Reference dataset provenance:
#   dataset : https://virtualcellmodels.cziscience.com/dataset/genome-scale-tcell-perturb-seq
#   code    : https://github.com/emdann/GWT_perturbseq_analysis_2025
#   preprint: https://www.biorxiv.org/content/10.64898/2025.12.23.696273v1
#   ~22M cells, 4 donors, conditions {Rest, Stim8hr, Stim48hr}, 26,504 sgRNAs / 12,654 genes + NTC.
#
# The pipeline consumes the PSEUDOBULK profile matrix (one profile per 10xrun x donor x condition x
# guide); the per-cell DE_stats objects are optional and large. Set DL_DE=1 to also fetch them.
set -euo pipefail
DATA_DIR="${DATA_DIR:-./data}"
OUT_DIR="$DATA_DIR/perturbseq"
mkdir -p "$OUT_DIR"

# --- reference CD4 T-cell dataset (override via env for other diseases) ---
PERTURB_URL="${PERTURB_URL:?set PERTURB_URL to the pseudobulk h5ad download URL (see dataset page)}"
PERTURB_OUT="${PERTURB_OUT:-$OUT_DIR/perturbseq.pseudobulk_merged.h5ad}"

echo "[perturbseq] $PERTURB_URL"
echo "         --> $PERTURB_OUT"
if [ ! -s "$PERTURB_OUT" ]; then
  curl -L --retry 5 --retry-delay 10 -C - -o "$PERTURB_OUT" "$PERTURB_URL"
else
  echo "         already present, skipping"
fi

if [ "${DL_DE:-0}" = "1" ]; then
  : "${DE_STATS_URL:?set DE_STATS_URL to fetch the optional per-cell DE_stats h5ad}"
  curl -L --retry 5 --retry-delay 10 -C - -o "$OUT_DIR/perturbseq.DE_stats.h5ad" "$DE_STATS_URL"
fi
echo "[perturbseq] done: $(du -h "$PERTURB_OUT" | cut -f1)"
