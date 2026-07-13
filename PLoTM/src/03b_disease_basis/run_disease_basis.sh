#!/bin/bash
# run_disease_basis.sh <config.yaml>  — OPTIONAL Stage-A replacement: learn cNMF programs from
# DISEASE single-cell data and project each Perturb-seq condition onto that fixed basis.
#
# Replaces the native `cnmf` phase for every condition CT. After it runs, the CT paths
#   <RESULTS>/cNMF/<CT>/test1/test1.{gene_spectra_score,gene_spectra_tpm,usages}.k_<Kd>.dt_0_4.*
# hold the DISEASE programs (identity, shared) + the condition's PROJECTED usages, so Stage B / C /
# Fig 4C / Fig 5a run unchanged with K = disease_singlecell.K.
#
# Called by drivers/run_disease.sh when the config has a disease_singlecell block. Standalone:
#   PY=... CNMF=... RSCRIPT=... RESULTS_ROOT=... DATA_DIR=... bash run_disease_basis.sh config/RA_singlecell.yaml
set -uo pipefail
CFG="${1:?usage: run_disease_basis.sh <config.yaml>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SRC="$HERE/src"; SC="$SRC/03b_disease_basis"
PY="${PY:?set PY to the cnmf-env python}"
CNMF="${CNMF:?set CNMF to the cnmf CLI}"
RSCRIPT="${RSCRIPT:?set RSCRIPT to an R env with Seurat (for .rds) }"
export DATA_DIR="${DATA_DIR:-$HERE/data}"
export RESULTS_ROOT="${RESULTS_ROOT:-$HERE/results}"
export PYTHONNOUSERSITE=1
cfg() { "$PY" "$SRC/cfg_get.py" "$CFG" "$1"; }

ID=$(cfg disease.id)
DCT=$(cfg disease_singlecell.ct_name); DCT="${DCT:-${ID}_sc}"
KD=$(cfg disease_singlecell.K)
NCORES=$(cfg compute.ncores)
RDS=$(cfg disease_singlecell.rds 2>/dev/null || true)
mapfile -t CONDS < <(cfg perturbseq.conditions)
ct_of(){ echo "${ID}_$(echo "$1" | tr 'A-Z' 'a-z')_pseudobulk"; }

echo "== disease-basis: DCT=$DCT K=$KD conds=${CONDS[*]} =="

# 0. .rds -> counts bundle (only if an rds path is configured and not already exported)
if [ -n "${RDS:-}" ] && [ "$RDS" != "None" ]; then
  LAYER=$(cfg disease_singlecell.counts_layer); LAYER="${LAYER:-counts}"
  PFX="$DATA_DIR/disease_sc/$DCT"
  if [ ! -s "$PFX.counts.mtx" ]; then
    echo "[sc/0] Seurat .rds -> counts bundle"
    "$RSCRIPT" "$SC/rds_to_counts.R" "$DATA_DIR/$RDS" "$LAYER" "$PFX" \
      || { RDS_ABS="$RDS"; "$RSCRIPT" "$SC/rds_to_counts.R" "$RDS_ABS" "$LAYER" "$PFX"; }
  fi
fi

# 1. build the disease cNMF input h5ad (intersect Perturb-measured genes + Stage-A QC)
echo "[sc/1] prep disease single-cell input"
"$PY" "$SC/prep_disease_singlecell.py" --config "$CFG" --data-dir "$DATA_DIR"

# 2. Stage A cNMF on the disease data -> the disease program basis (K=$KD)
echo "[sc/2] cNMF on disease data (K=$KD)"
NCORES="$NCORES" bash "$SRC/03_cnmf_stageA/run_cnmf.sh" "$DCT" "$KD" "$NCORES"

# 3. filter each Perturb-seq condition (Stage-A QC) then project onto the fixed disease basis
for c in "${CONDS[@]}"; do CT=$(ct_of "$c")
  FILT="$DATA_DIR/Perturbseq/filtered_data/$CT.h5ad"
  if [ ! -s "$FILT" ]; then
    echo "[sc/3a] filter_cells $CT"
    ( cd "$HERE" && "$PY" "$SRC/03_cnmf_stageA/filter_cells.py" "$CT" )
  fi
  echo "[sc/3b] project $CT onto $DCT"
  "$PY" "$SC/project_onto_disease.py" \
      --disease-ct "$DCT" --disease-k "$KD" \
      --condition-h5ad "$FILT" \
      --ct "$CT" --metadata "$DATA_DIR/Perturbseq/metadata/${CT}_metadata.csv" \
      --results "$RESULTS_ROOT"
done
echo "== disease-basis done: programs=$DCT (K=$KD); usages projected for ${#CONDS[@]} conditions =="
