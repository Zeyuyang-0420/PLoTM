#!/bin/bash
# run_disease.sh <config.yaml> <phase>   — top-level, config-driven orchestrator.
#
# Phases (each depends on the previous; use `all` for the lot):
#   download   perturb-seq h5ad + LoF burden summary-stats from both sources
#   posterior  GeneBayes gene-level LoF gamma per trait (GPU)  [+ QC gates]
#   input      build the cNMF input h5ad + covariate metadata per condition
#   cnmf       Stage A: cNMF programs per condition
#   stageB     Stage B: regulatory effect sizes beta_x(P) per condition
#   thresholds per-trait |gamma| cutoff (q99) + QC flags
#   burden     Stage C: program/regulator burden per (trait, condition)  -> Fig 4C inputs
#   figures    Fig 4C + Fig 5a (obs + permutation) per (trait, condition)
#   annotate   program annotation (GO/hallmark/markers) per condition
#   analysis   condition comparison: tidy table, meta-analysis, networks, priority checks
#
# The heavy per-(trait,condition) loop (burden/figures) is delegated to run_batch.sh; this script
# wires the UPSTREAM disease-specific stages and reads all knobs from the YAML config.
#
# This is the interface an AGENT drives for a new disease: point the config at a new perturb-seq
# dataset + new phenotype list and call `all`. See docs/AGENT_GUIDE.md.
#
# INTEGRATION STATUS (be honest with yourself, agent — do not overclaim automation):
#   CONFIG-DRIVEN FROM SCRATCH : input · cnmf · stageB · thresholds
#       (make_pseudobulk_input.py / run_cnmf.sh / regulatory_effectsize.R / compute_thresholds.py all
#        take their inputs from this config or from env vars this driver sets.)
#   ENV-OVERRIDABLE REFERENCE  : burden · figures · analysis
#       (run_batch.sh + src/08_condition_analysis/* default to the RA reference paths; export
#        RESULTS_ROOT/DATA_ROOT and edit the trait_groups for a new disease. They RUN as-is for RA.)
#   NEEDS PER-DISEASE IDENTIFIERS : download + posterior
#       (fetch_*_lof.py ship the RA phenotype panel; burden summary stats cannot be auto-discovered
#        from a phenotype name — supply each GeneBass phenocode / Backman GWAS-Catalog accession.)
set -uo pipefail
CFG="${1:?usage: run_disease.sh <config.yaml> <phase>}"
PHASE="${2:?usage: run_disease.sh <config.yaml> <phase>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="$HERE/src"
PY="${PY:?set PY to the analysis python (pyyaml, scanpy, pandas)}"
RSCRIPT="${RSCRIPT:?set RSCRIPT to the R env}"

# --- read config into shell vars via a tiny python helper ---
cfg() { "$PY" "$SRC/cfg_get.py" "$CFG" "$1"; }
ID=$(cfg disease.id)
mapfile -t CONDS < <(cfg perturbseq.conditions)
mapfile -t TRAITS < <(cfg lof.traits.names)
K=$(cfg params.K)
# OPTIONAL disease single-cell cNMF basis: if config has a disease_singlecell block, programs are
# learned from disease single-cell data and Perturb-seq conditions are PROJECTED onto that fixed
# basis. Downstream K becomes the disease K. Absent -> pure Perturb-seq (native cNMF), K=params.K.
DSC_CT=$(cfg disease_singlecell.ct_name 2>/dev/null || true)
if [ -n "${DSC_CT:-}" ] && [ "$DSC_CT" != "None" ]; then
  DISEASE_BASIS=1; K=$(cfg disease_singlecell.K)
  echo "== MODE: disease single-cell basis ($DSC_CT, K=$K) — Perturb conditions will be projected =="
else
  DISEASE_BASIS=0
  echo "== MODE: pure Perturb-seq (native cNMF, K=$K) =="
fi
export DATA_DIR="${DATA_DIR:-$HERE/data}"
export RESULTS_ROOT="${RESULTS_ROOT:-$HERE/results/$ID}"
export LOF_THRESH_MODE=$(cfg params.lof_threshold)
NCORES=$(cfg compute.ncores); THROTTLE=$(cfg compute.throttle)
ct_of(){ echo "${ID}_$(echo "$1" | tr 'A-Z' 'a-z')_pseudobulk"; }

echo "== run_disease id=$ID phase=$PHASE conds=${CONDS[*]} traits=${#TRAITS[@]} K=$K =="
has(){ [[ "$PHASE" == "all" || "$PHASE" == "$1" ]]; }

if has download; then
  echo "[download] perturb-seq + LoF (genebass + backman)"
  PERTURB_URL=$(cfg perturbseq.url) bash "$SRC/00_download/download_perturbseq.sh"
  "$PY" "$SRC/00_download/fetch_genebass_lof.py" --config "$CFG"
  "$PY" "$SRC/00_download/fetch_backman_lof.py"  --config "$CFG"
  "$PY" "$SRC/00_download/prepare_lof_inputs_genebass.py" --config "$CFG"
  "$PY" "$SRC/00_download/prepare_lof_inputs_backman.py"  --config "$CFG"
fi

if has posterior; then
  echo "[posterior] GeneBayes gamma per trait (GPU)"
  bash "$SRC/01_lof_posterior/run_genebayes_converged.sh" "${TRAITS[@]}"
  "$PY" "$SRC/01_lof_posterior/qc_posteriors.py" --config "$CFG"
fi

if has input; then
  PB=$(cfg perturbseq.pseudobulk_h5ad)
  for c in "${CONDS[@]}"; do
    echo "[input] $c -> $(ct_of "$c")"
    "$PY" "$SRC/02_input_prep/make_pseudobulk_input.py" \
        --input "$PB" --condition "$c" --ct "$(ct_of "$c")"
  done
fi

if has cnmf; then
  if [ "$DISEASE_BASIS" = "1" ]; then
    # learn programs from disease single-cell data + project every Perturb condition onto it
    NCORES="$NCORES" bash "$SRC/03b_disease_basis/run_disease_basis.sh" "$CFG"
  else
    for c in "${CONDS[@]}"; do CT=$(ct_of "$c")
      echo "[cnmf] Stage A: $CT"
      "$PY" "$SRC/03_cnmf_stageA/filter_cells.py" "$CT"
      NCORES="$NCORES" bash "$SRC/03_cnmf_stageA/run_cnmf.sh" "$CT" "$K" "$NCORES"
    done
  fi
fi

if has stageB; then
  MINP=$(cfg params.stageB_min_profiles); MINP="${MINP:-4}"
  for c in "${CONDS[@]}"; do CT=$(ct_of "$c")
    echo "[stageB] regulatory effects: $CT (K=$K programs)"
    # pass K (and min-profiles) so Stage B matches the disease K; non-estimable programs are guarded.
    for i in $(seq 1 "$K"); do "$RSCRIPT" "$SRC/04_regulatory_stageB/regulatory_effectsize.R" "$CT" "$i" "$K" "$MINP" & done; wait
  done
fi

if has thresholds; then
  echo "[thresholds] per-trait q99 + QC flags"
  POSTERIOR_DIR="$DATA_DIR/lof/GeneBayes_posterior" THRESH_OUT="$RESULTS_ROOT/lof_thresholds.tsv" \
    "$PY" "$HERE/drivers/compute_thresholds.py"
fi

# ---- burden + figures: delegate the (trait x condition) matrix to run_batch.sh ----
if has burden || has figures; then
  P=(); has burden && P+=(C fig4c); has figures && P+=(fig5obs fig5perm fig5a)
  echo "[batch] phases: ${P[*]}"
  NPROC="$THROTTLE" STAGES="${P[*]}" bash "$HERE/drivers/run_batch.sh" "${P[@]}"
fi

if has annotate; then
  for c in "${CONDS[@]}"; do CT=$(ct_of "$c")
    echo "[annotate] $CT"
    "$RSCRIPT" "$SRC/06_annotation/annotate_programs.R" "$CT" || echo "  (annotation optional; needs clusterProfiler)"
  done
fi

if has analysis; then
  echo "[analysis] condition comparison + meta-analysis + networks"
  cd "$SRC/08_condition_analysis"
  for s in build_tidy meta_analysis condition_interaction priority_checks overlap_permutation make_figures; do
    "$PY" "$s.py" --config "$CFG" || echo "  $s.py needs the tidy table / results present"
  done
fi
echo "== run_disease $PHASE done -> $RESULTS_ROOT =="
