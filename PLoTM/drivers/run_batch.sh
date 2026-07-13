#!/bin/bash
# run_batch.sh <PHASE> [TRAIT ...]
#
# Batch driver: 8 converged GeneBayes LoF posteriors x 2 conditions -> Fig 4C + Fig 5a,
# all against the SAME cNMF program matrix + Stage-B betas (which are LoF-independent).
#
# Phases (run in this order):
#   prep      symlink posterior into DATA_ROOT, build the symlinked data/ work tree, validate
#   C         Stage-C burden test (60 programs x 2 conditions per trait), global process pool
#   fig4c     Figure 4C per (trait, condition)
#   fig5obs   Fig-5 observed pass  -- LOW concurrency on purpose (regsubsets is memory-bandwidth bound)
#   fig5perm  Fig-5 permutation null, global process pool of shards
#   fig5a     Figure 5a assembly + permutation P
#   collect   gather figures/tables into pipeline_results/<trait>/
#
# Concurrency: NPROC (default 80) leaves headroom for the GeneBayes GPU job still running.
# Per-trait |gamma| threshold is read from lof_thresholds.tsv (99th pct of |gamma|), NOT a
# fixed 0.03: Backman gamma is ~60x the Genebass scale, so a fixed cutoff would select
# every gene for Backman and leave an empty background.

set -uo pipefail
PHASE="${1:?usage: run_batch.sh <prep|C|fig4c|fig5obs|fig5perm|fig5a|collect> [TRAIT...]}"
shift

PIPE=/mnt/scratch/ZY2/Hackathon/RA_cNMF_LoF_pipeline
ROOT=/mnt/scratch/ZY2/Hackathon/pipeline_results
POST=/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/RA_LoF/data/LoF/GeneBayes_posterior
export DATA_ROOT=/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data
export RESULTS_ROOT=/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis
RSCRIPT=/mnt/scratch/ZY2/.envs/env_r_core/bin/Rscript
CTS=(GWCD4i_rest_pseudobulk GWCD4i_stim48_pseudobulk)
K=60
NPROC="${NPROC:-80}"
PN="${PN:-5}"; RN="${RN:-3}"; TOP="${TOP:-200}"
NPERM="${NPERM:-5000}"; NPERM_OBS="${NPERM_OBS:-10000}"; NSHARD="${NSHARD:-10}"

# traits: args, else every READY row in lof_thresholds.tsv
if [ "$#" -gt 0 ]; then TRAITS=("$@"); else
  mapfile -t TRAITS < <(awk -F'\t' 'NR>1 && $2=="READY"{print $1}' "$ROOT/lof_thresholds.tsv")
fi
thresh_of(){ awk -F'\t' -v t="$1" 'NR>1 && $1==t{print $4}' "$ROOT/lof_thresholds.tsv"; }
lab_of(){ [[ "$1" == *rest* ]] && echo "GWCD4i Rest" || echo "GWCD4i Stim48hr"; }

echo "== phase=$PHASE traits=${#TRAITS[@]} nproc=$NPROC =="
[ "${#TRAITS[@]}" -eq 0 ] && { echo "no traits"; exit 1; }

case "$PHASE" in

prep)
  for t in "${TRAITS[@]}"; do
    src="$POST/${t}.per_gene_estimates.tsv"
    [ -s "$src" ] || { echo "  SKIP $t (posterior not ready)"; continue; }
    ln -sfn "$src" "$DATA_ROOT/LoF/GeneBayes_posterior/${t}.per_gene_estimates.tsv"
    th=$(thresh_of "$t")
    for ct in "${CTS[@]}"; do
      WORKDIR="$ROOT/$t/$ct" STAGES=none LOF_THRESH="$th" \
        bash "$PIPE/run_trait.sh" "$ct" "${t}.per_gene_estimates.tsv" "$K" >/dev/null || { echo "  FAIL prep $t/$ct"; exit 1; }
    done
    echo "  prep $t  (thresh=$th)"
  done
  ;;

C)
  TASKS="$ROOT/logs/tasks.C.$$"; : > "$TASKS"
  for t in "${TRAITS[@]}"; do for ct in "${CTS[@]}"; do
    wd="$ROOT/$t/$ct"; mkdir -p "$wd/data/Perturbseq/trait_association/$ct/ProgramLevel/tmp_parallel" "$wd/logs"
    for i in $(seq 1 $K); do echo -e "$wd\t$ct\t${t}.per_gene_estimates.tsv\t$i" >> "$TASKS"; done
  done; done
  echo "  $(wc -l < "$TASKS") program-tasks -> pool of $NPROC (skips programs already done)"
  # RESUMABLE: skip any program whose per-program output already exists and is non-empty.
  # The box is shared (a separate 80-proc fig5 job + the GeneBayes GPU driver), so this phase
  # is run in waves at whatever concurrency is free; re-invoking must not redo finished work.
  < "$TASKS" xargs -P "$NPROC" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r wd ct lof i <<< "{}"
    out="$wd/data/Perturbseq/trait_association/$ct/ProgramLevel/tmp_parallel/programs_enrichment_K'"$K"'_${lof}.p$i"
    reg="$wd/data/Perturbseq/trait_association/$ct/ProgramLevel/tmp_parallel/regulators_enrichment_K'"$K"'_${lof}.p$i"
    # complete == header + >=1 data row in BOTH files. A wave killed mid-flight can leave a
    # truncated file; -s alone would then wrongly mark it done.
    if [ -f "$out" ] && [ -f "$reg" ] \
       && [ "$(wc -l < "$out")" -ge 2 ] && [ "$(wc -l < "$reg")" -ge 2 ]; then exit 0; fi
    cd "$wd" && '"$RSCRIPT"' '"$PIPE"'/scripts/5_burden_pseudobulk_oneprogram.R "$lof" '"$K"' "$ct" "$i" \
      > "logs/burden.p$i.log" 2>&1 || echo "FAIL $wd p$i"'
  rm -f "$TASKS"
  # concatenate per-program shards into the standard enrichment tables
  for t in "${TRAITS[@]}"; do for ct in "${CTS[@]}"; do
    PL="$ROOT/$t/$ct/data/Perturbseq/trait_association/$ct/ProgramLevel"
    for kind in regulators programs; do
      f="$PL/${kind}_enrichment_K${K}_${t}.per_gene_estimates.tsv"
      head -1 "$PL/tmp_parallel/${kind}_enrichment_K${K}_${t}.per_gene_estimates.tsv.p1" > "$f" 2>/dev/null || continue
      for i in $(seq 1 $K); do tail -n +2 "$PL/tmp_parallel/${kind}_enrichment_K${K}_${t}.per_gene_estimates.tsv.p$i" >> "$f" 2>/dev/null; done
    done
    n=$(( $(wc -l < "$PL/programs_enrichment_K${K}_${t}.per_gene_estimates.tsv" 2>/dev/null || echo 1) - 1 ))
    echo "  $t/$ct: $n/60 programs"
  done; done
  ;;

fig4c)
  for t in "${TRAITS[@]}"; do for ct in "${CTS[@]}"; do
    wd="$ROOT/$t/$ct"; PL="data/Perturbseq/trait_association/$ct/ProgramLevel"
    # LABEL is the CONDITION only: plot_fig4_burden.R already appends the trait to both the
    # title ("Fig 4C | <LABEL> | <trait>") and the filename.
    ( cd "$wd" && "$RSCRIPT" "$PIPE/scripts/plot_fig4_burden.R" "${t}.per_gene_estimates.tsv" "$K" "$PL" figures "$(lab_of "$ct")" ) \
      2>&1 | grep -E "Bonferroni-significant|saved" | sed "s|^|  $t/$ct: |"
  done; done
  ;;

fig5obs)
  TASKS="$ROOT/logs/tasks.obs.$$"; : > "$TASKS"
  for t in "${TRAITS[@]}"; do th=$(thresh_of "$t"); for ct in "${CTS[@]}"; do
    echo -e "$ROOT/$t/$ct\t$ct\t${t}.per_gene_estimates.tsv\t$th" >> "$TASKS"; done; done
  echo "  $(wc -l < "$TASKS") observed passes (concurrency $(wc -l < "$TASKS"), each 1 proc)"
  < "$TASKS" xargs -P 16 -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r wd ct lof th <<< "{}"
    cd "$wd" && '"$RSCRIPT"' '"$PIPE"'/scripts/fig5_permutation_step1_pseudobulk.R "$ct" "$lof" '"$PN $RN $TOP"' "$th" obs 1 0 '"$NPERM_OBS"' fig5_out \
      > logs/fig5.obs.log 2>&1 || echo "FAIL obs $wd"'
  rm -f "$TASKS"
  for t in "${TRAITS[@]}"; do for ct in "${CTS[@]}"; do
    grep -h "OBSERVED programs=" "$ROOT/$t/$ct/logs/fig5.obs.log" 2>/dev/null | tail -1 | sed "s|^|  $t/$ct: |" \
      || echo "  $t/$ct: NO OBSERVED LINE"; done; done
  ;;

fig5perm)
  TASKS="$ROOT/logs/tasks.perm.$$"; : > "$TASKS"
  per=$(( (NPERM + NSHARD - 1) / NSHARD ))
  for t in "${TRAITS[@]}"; do th=$(thresh_of "$t"); for ct in "${CTS[@]}"; do
    s=1; while [ "$s" -le "$NPERM" ]; do e=$(( s + per - 1 )); [ "$e" -gt "$NPERM" ] && e="$NPERM"
      echo -e "$ROOT/$t/$ct\t$ct\t${t}.per_gene_estimates.tsv\t$th\t$s\t$e" >> "$TASKS"; s=$(( e + 1 )); done
  done; done
  echo "  $(wc -l < "$TASKS") shards x $per perms -> pool of $NPROC (skips completed shards)"
  # RESUMABLE: a shard is done iff its log records "DONE (perm".
  < "$TASKS" xargs -P "$NPROC" -I{} bash -c '
    IFS=$'"'"'\t'"'"' read -r wd ct lof th s e <<< "{}"
    lg="$wd/logs/fig5.perm_${s}_${e}.log"
    grep -q "DONE (perm" "$lg" 2>/dev/null && exit 0
    cd "$wd" && '"$RSCRIPT"' '"$PIPE"'/scripts/fig5_permutation_step1_pseudobulk.R "$ct" "$lof" '"$PN $RN $TOP"' "$th" perm "$s" "$e" '"$NPERM_OBS"' fig5_out \
      > "$lg" 2>&1 || echo "FAIL perm $wd $s-$e"'
  rm -f "$TASKS"
  for t in "${TRAITS[@]}"; do for ct in "${CTS[@]}"; do
    d=$(grep -l 'DONE (perm' "$ROOT/$t/$ct/logs/"fig5.perm_*.log 2>/dev/null | wc -l); echo "  $t/$ct: $d/$NSHARD shards"; done; done
  ;;

fig5a)
  for t in "${TRAITS[@]}"; do th=$(thresh_of "$t"); for ct in "${CTS[@]}"; do
    ( cd "$ROOT/$t/$ct" && "$RSCRIPT" "$PIPE/scripts/fig5_step2_assemble.R" "$ct" "${t}.per_gene_estimates.tsv" \
        "$PN" "$RN" "$TOP" "$th" fig5_out figures 2000 ) 2>&1 | grep -E "^\[step2\]" | sed "s|^|  $t/$ct: |"
  done; done
  ;;

collect)
  for t in "${TRAITS[@]}"; do
    out="$ROOT/$t"; mkdir -p "$out/figures" "$out/tables"
    for ct in "${CTS[@]}"; do
      short=$([[ "$ct" == *rest* ]] && echo rest || echo stim48hr)
      for f in "$ROOT/$t/$ct/figures/"*.png "$ROOT/$t/$ct/figures/"*.pdf; do
        [ -e "$f" ] || continue; cp -f "$f" "$out/figures/${short}__$(basename "$f")"; done
      PL="$ROOT/$t/$ct/data/Perturbseq/trait_association/$ct/ProgramLevel"
      for f in "$PL/programs_enrichment_K${K}_${t}.per_gene_estimates.tsv" "$PL/regulators_enrichment_K${K}_${t}.per_gene_estimates.tsv"; do
        [ -e "$f" ] && cp -f "$f" "$out/tables/${short}__$(basename "$f")"; done
      for f in "$ROOT/$t/$ct/fig5_out/$ct/P${TOP}/"*.txt "$ROOT/$t/$ct/fig5_out/$ct/P${TOP}/"*.tsv; do
        [ -e "$f" ] || continue; cp -f "$f" "$out/tables/${short}__$(basename "$f")"; done
    done
    echo "  $t: $(ls "$out/figures" 2>/dev/null | wc -l) figures, $(ls "$out/tables" 2>/dev/null | wc -l) tables"
  done
  ;;

*) echo "unknown phase $PHASE"; exit 1;;
esac
echo "== $PHASE done =="
