#!/bin/bash
# run_trait.sh <CT> <LOF_FILE> [K=60]
#
# End-to-end: given ONE already-computed cNMF matrix (+ Stage-B betas) and ONE LoF trait, produce
#   Figure 4C  (program-burden vs regulator-burden scatter)
#   Figure 5a  (regulators -> programs -> trait map + permutation P)
#
# The cNMF matrix and the Stage-B betas are LoF-INDEPENDENT: compute them once (see prerequisites/)
# and re-run this script for each trait. The `data/` tree is assembled from SYMLINKS so every trait
# reads the identical cNMF matrix.
#
#   CT        e.g. GWCD4i_rest_pseudobulk
#   LOF_FILE  e.g. Genebass_RA_M06.per_gene_estimates.tsv   (must sit in $DATA_ROOT/LoF/GeneBayes_posterior/)
#
# Required env (or edit the defaults):
#   RESULTS_ROOT  dir containing  cNMF/<CT>/test1/  and  cNMF_regulation/<CT>/
#   DATA_ROOT     dir containing  gencode_v41_gname_gid_ALL_sorted_onlyID, shet_10bins.txt,
#                                 LoF/GeneBayes_posterior/<LOF_FILE>
#   WORKDIR       scratch dir for the symlinked data/ tree + outputs   (default: ./work)
#   RSCRIPT       Rscript with: data.table, dplyr, leaps, ggplot2, ggrepel
#
# Tunables: PN RN TOP LOF_THRESH NPERM NPERM_OBS NSHARD THROTTLE  (see README)
# Flags:    STAGES="C fig4c fig5obs fig5perm fig5a"   to run a subset
#
# NOTE (learned the hard way):
#   * Run the Fig-5 OBSERVED pass ALONE, before the permutation shards. With ~100 concurrent R procs
#     `leaps::regsubsets` becomes memory-bandwidth bound (~33 s/perm vs ~9 s), and the observed pass
#     — which gates step 2 — gets starved. This script enforces that ordering.
#   * Never `pkill -f` on a pattern that also matches this wrapper's own command line.
#   * Do not use foreground `sleep` in the harness.

set -uo pipefail
CT="${1:?usage: run_trait.sh <CT> <LOF_FILE> [K]}"
LOF_FILE="${2:?usage: run_trait.sh <CT> <LOF_FILE> [K]}"
K="${3:-60}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_ROOT="${RESULTS_ROOT:-/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis}"
DATA_ROOT="${DATA_ROOT:-/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data}"
WORKDIR="${WORKDIR:-$HERE/work}"
RSCRIPT="${RSCRIPT:-/mnt/scratch/ZY2/.envs/env_r_core/bin/Rscript}"

PN="${PN:-5}"; RN="${RN:-3}"; TOP="${TOP:-200}"; LOF_THRESH="${LOF_THRESH:-0.03}"
NPERM="${NPERM:-5000}"; NPERM_OBS="${NPERM_OBS:-10000}"; NSHARD="${NSHARD:-26}"; THROTTLE="${THROTTLE:-40}"
STAGES="${STAGES:-C fig4c fig5obs fig5perm fig5a}"
TRAIT="${LOF_FILE%.per_gene_estimates.tsv}"
LABEL="$( [[ "$CT" == *rest* ]] && echo "GWCD4i Rest" || echo "GWCD4i Stim48hr" )"
has(){ [[ " $STAGES " == *" $1 "* ]]; }

echo "== run_trait CT=$CT TRAIT=$TRAIT K=$K PN=$PN RN=$RN TOP=$TOP LOF=$LOF_THRESH NPERM=$NPERM =="

## ---------- 0. validate inputs, assemble the symlinked data/ tree ----------
need(){ [[ -e "$1" ]] || { echo "MISSING: $1"; exit 1; }; }
need "$RESULTS_ROOT/cNMF/$CT/test1/test1.gene_spectra_score.k_${K}.dt_0_4.txt"
need "$RESULTS_ROOT/cNMF_regulation/$CT/K${K}_program1_perturb_effects.txt"
need "$DATA_ROOT/gencode_v41_gname_gid_ALL_sorted_onlyID"
need "$DATA_ROOT/shet_10bins.txt"
need "$DATA_ROOT/LoF/GeneBayes_posterior/$LOF_FILE"

mkdir -p "$WORKDIR/data/Perturbseq/cNMF" "$WORKDIR/data/Perturbseq/cNMF_regulation" \
         "$WORKDIR/data/Perturbseq/trait_association/$CT/ProgramLevel" "$WORKDIR/figures" "$WORKDIR/logs"
ln -sfn "$RESULTS_ROOT/cNMF/$CT"            "$WORKDIR/data/Perturbseq/cNMF/$CT"
ln -sfn "$RESULTS_ROOT/cNMF_regulation/$CT" "$WORKDIR/data/Perturbseq/cNMF_regulation/$CT"
ln -sfn "$DATA_ROOT/LoF"                    "$WORKDIR/data/LoF"
ln -sf  "$DATA_ROOT/gencode_v41_gname_gid_ALL_sorted_onlyID" "$WORKDIR/data/"
ln -sf  "$DATA_ROOT/shet_10bins.txt"                         "$WORKDIR/data/"
cd "$WORKDIR"
export RESULTS_ROOT DATA_ROOT
PL="data/Perturbseq/trait_association/$CT/ProgramLevel"

## ---------- 1. Stage C: burden test (Fig-4C inputs) ----------
if has C; then
  echo "[1] Stage C burden: $K programs in parallel (throttle $THROTTLE)"
  mkdir -p "$PL/tmp_parallel"; pids=()
  for i in $(seq 1 "$K"); do
    "$RSCRIPT" "$HERE/scripts/5_burden_pseudobulk_oneprogram.R" "$LOF_FILE" "$K" "$CT" "$i" \
        > "logs/burden.$TRAIT.p$i.log" 2>&1 &
    pids+=($!); (( ${#pids[@]} >= THROTTLE )) && { wait "${pids[0]}"; pids=("${pids[@]:1}"); }
  done
  wait
  for kind in regulators programs; do
    f="$PL/${kind}_enrichment_K${K}_${LOF_FILE}"
    head -1 "$PL/tmp_parallel/${kind}_enrichment_K${K}_${LOF_FILE}.p1" > "$f"
    for i in $(seq 1 "$K"); do tail -n +2 "$PL/tmp_parallel/${kind}_enrichment_K${K}_${LOF_FILE}.p$i" >> "$f" 2>/dev/null; done
    echo "    $(basename "$f"): $(( $(wc -l < "$f") - 1 )) programs"
  done
fi

## ---------- 2. Figure 4C ----------
if has fig4c; then
  echo "[2] Figure 4C"
  "$RSCRIPT" "$HERE/scripts/plot_fig4_burden.R" "$LOF_FILE" "$K" "$PL" "figures" "$LABEL" 2>&1 \
    | grep -E "Bonferroni-significant|labelled|saved" || true
fi

## ---------- 3. Fig-5 observed pass (MUST run alone: see header note) ----------
if has fig5obs; then
  echo "[3] Fig5 observed pass (NPERM_OBS=$NPERM_OBS) — running alone"
  "$RSCRIPT" "$HERE/scripts/fig5_permutation_step1_pseudobulk.R" "$CT" "$LOF_FILE" \
      "$PN" "$RN" "$TOP" "$LOF_THRESH" obs 1 0 "$NPERM_OBS" fig5_out 2>&1 | grep -E "^\[fig5\]" || true
fi

## ---------- 4. Fig-5 permutation null (sharded; equivalent to the serial loop) ----------
if has fig5perm; then
  echo "[4] Fig5 permutation null: $NPERM perms in $NSHARD shards"
  per=$(( (NPERM + NSHARD - 1) / NSHARD )); s=1; pids=()
  while [ "$s" -le "$NPERM" ]; do
    e=$(( s + per - 1 )); [ "$e" -gt "$NPERM" ] && e="$NPERM"
    "$RSCRIPT" "$HERE/scripts/fig5_permutation_step1_pseudobulk.R" "$CT" "$LOF_FILE" \
        "$PN" "$RN" "$TOP" "$LOF_THRESH" perm "$s" "$e" "$NPERM_OBS" fig5_out \
        > "logs/fig5.$TRAIT.perm_${s}_${e}.log" 2>&1 &
    pids+=($!); s=$(( e + 1 ))
  done
  wait
  echo "    shards done: $(grep -l 'DONE (perm' logs/fig5.$TRAIT.perm_*.log 2>/dev/null | wc -l)/$NSHARD"
fi

## ---------- 5. Figure 5a (+ permutation histogram, concordance scatter) ----------
if has fig5a; then
  echo "[5] Figure 5a assembly"
  "$RSCRIPT" "$HERE/scripts/fig5_step2_assemble.R" "$CT" "$LOF_FILE" "$PN" "$RN" "$TOP" "$LOF_THRESH" \
      fig5_out figures 2000 2>&1 | grep -E "^\[step2\]" || true
fi

echo "== done. outputs in $WORKDIR/figures and $WORKDIR/fig5_out/$CT/P$TOP =="
