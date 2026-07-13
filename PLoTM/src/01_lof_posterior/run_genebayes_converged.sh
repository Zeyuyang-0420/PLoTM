#!/bin/bash
# 10_run_converged.sh [TRAIT...]
# Produce the *fully paper-faithful* posterior: --total_iterations 1000 (the paper's value).
# NGBoost early stopping (patience 10) fires around iteration ~140-150 and the script feeds
# best_val_loss_itr to the posterior, so this costs ~70 min/trait, not hours.
# Outputs suffixed _conv so earlier runs are preserved.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /mnt/scratch/ZY2/.envs/env_genebayes
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=0

ROOT="/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/RA_LoF"
SCRIPT="/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/source_code/data_preparation/GeneBayes/genebayes.lof.gene_level.py"
OUTDIR="$ROOT/genebayes_run"; POSTDIR="$ROOT/data/LoF/GeneBayes_posterior"
mkdir -p "$OUTDIR" "$POSTDIR"

for NAME in "$@"; do
  TAG="${NAME}_conv"
  if [ -s "$POSTDIR/${TAG}.per_gene_estimates.tsv" ]; then echo "SKIP $TAG"; continue; fi
  echo "======== [$(date +%H:%M:%S)] GeneBayes 1000-iter (early-stops): $NAME ========"
  python "$SCRIPT" \
    --response    "$ROOT/traits/${NAME}.summary_statistics.csv" \
    --features    "$ROOT/features/gene_features_ensg_numeric.tsv" \
    --train_genes "$ROOT/genelists/${NAME}.all.txt" \
    --val_genes   "$ROOT/genelists/${NAME}.all.txt" \
    --out         "$OUTDIR/${TAG}" \
    --integration_lb -3.99 --integration_ub 3.99 \
    --batch_size 100 --n_integration_pts 50000 \
    --lr 0.3 --total_iterations 1000 \
    --n_trees_per_iteration 2 --reg_alpha 2 --reg_lambda 2 \
    --subsample 0.8 --min_child_weight 3 --max_depth 3 \
    > "$OUTDIR/${TAG}.run.log" 2>&1
  cp "$OUTDIR/${TAG}.per_gene_estimates.tsv" "$POSTDIR/${TAG}.per_gene_estimates.tsv"
  echo "CONV_TRAIT_DONE $TAG [$(date +%H:%M:%S)]"
done
echo "CONV_ALL_DONE"
