#!/bin/bash
# run_cnmf.sh — Stage A: cNMF consensus program discovery (Kotliar et al.), paper-faithful settings.
#
#   run_cnmf.sh <CT> [K=60] [NCORES]
#
# Verbatim cNMF CLI: -k K --n-iter 100 --seed 14 --numgenes 2000, consensus at
# --local-density-threshold 0.4  (written as dt_0_4). The only harness change vs the paper's SLURM
# scripts is that the factorize array is a local loop over worker indices.
#
# Inputs : data/perturbseq/filtered_data/<CT>.h5ad  (from filter_cells.py)
# Outputs: $RESULTS_ROOT/cNMF/<CT>/test1/test1.gene_spectra_score.k_<K>.dt_0_4.txt   (60 x genes)
#          $RESULTS_ROOT/cNMF/<CT>/test1/test1.usages.k_<K>.dt_0_4.consensus.txt      (profiles x 60)
set -euo pipefail
CT="${1:?usage: run_cnmf.sh <CT> [K] [NCORES]}"
K="${2:-60}"
NCORES="${3:-${NCORES:-64}}"
PY="${PY:?set PY to the cnmf-env python}"
CNMF="${CNMF:?set CNMF to the cnmf CLI}"
DATA_DIR="${DATA_DIR:-./data}"
RESULTS_ROOT="${RESULTS_ROOT:-./results}"

# thread pinning keeps peak RAM bounded during the 100-replicate factorization
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONNOUSERSITE=1

OUT="$RESULTS_ROOT/cNMF/$CT"
IN="$DATA_DIR/perturbseq/filtered_data/$CT.h5ad"
[ -s "$IN" ] || { echo "MISSING filtered input $IN (run filter_cells.py first)"; exit 1; }
mkdir -p "$OUT"

echo "[A/1] prepare (K=$K, 100 iters, seed 14, 2000 HVGs)"
"$CNMF" prepare --output-dir "$OUT" --name test1 -c "$IN" \
    -k "$K" --n-iter 100 --seed 14 --numgenes 2000 --total-workers "$NCORES"

echo "[A/2] factorize ($NCORES workers)"
for i in $(seq 0 $((NCORES-1))); do
    "$CNMF" factorize --output-dir "$OUT" --name test1 --worker-index "$i" --total-workers "$NCORES" &
done
wait

echo "[A/3] combine + consensus (dt=0.4)"
"$CNMF" combine --output-dir "$OUT" --name test1
"$CNMF" consensus --output-dir "$OUT" --name test1 --components "$K" \
    --local-density-threshold 0.4 --show-clustering
echo "== Stage A done: $(ls -1 "$OUT"/test1/*.consensus.txt 2>/dev/null | wc -l) consensus files -> $OUT/test1 =="
