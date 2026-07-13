#!/usr/bin/env python3
"""
13_ensemble_gamma.py — average γ over independent CONVERGED GeneBayes runs.

Why: `genebayes.lof.gene_level.py` sets no random seed and is stochastic (randperm batching,
XGBoost subsample=0.8, a 10k-sample MC prior_mean). Two independent converged runs of the
same trait on the same data agree at only r ≈ 0.70 (56/100 top-γ genes shared). γ is a noisy
but apparently unbiased estimator, so averaging k runs raises reliability per Spearman–Brown:
    rel(k) = k·r / (1 + (k−1)·r)      r=0.70:  k=2 → 0.82, k=3 → 0.87, k=5 → 0.92
Verified: the 2-run ensemble cross-trait γ correlation was predicted at +0.6161 and observed
at +0.6157 (error 4e-4).

Only *converged* runs are averaged (`_it300`, `_conv`). The 50-iteration posteriors are
under-converged (biased, not merely noisy) and are deliberately excluded.

Writes data/LoF/GeneBayes_posterior/<stem>_ens.per_gene_estimates.tsv with the columns the
pipeline reads (`ensg`, `post_mean`) plus `prior_mean` and `n_runs` for provenance.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
POST = ROOT / "data" / "LoF" / "GeneBayes_posterior"
CONVERGED_SUFFIXES = ["_it300", "_conv"]   # NOT the bare 50-iter run

STEMS = [
    "Genebass_RF", "Genebass_RA_custom", "Genebass_RA_M06",
    "Backman_2021_RA", "Backman_2021_RA_M06", "Backman_2021_RA_alt",
    "Backman_2021_RA_M05_seropos", "Backman_2021_RA_M060_seroneg",
]


def load(p):
    g, pr = {}, {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["post_mean"] in ("", "nan"):
            continue
        g[r["ensg"]] = float(r["post_mean"])
        pr[r["ensg"]] = float(r["prior_mean"])
    return g, pr


def sb(r, k):
    return k * r / (1 + (k - 1) * r)


for stem in STEMS:
    runs = [(POST / f"{stem}{s}.per_gene_estimates.tsv") for s in CONVERGED_SUFFIXES]
    runs = [p for p in runs if p.exists()]
    if not runs:
        print(f"[{stem}] no converged run found — skipping")
        continue
    gs, ps = zip(*(load(p) for p in runs))
    common = set(gs[0])
    for g in gs[1:]:
        common &= set(g)
    common = sorted(common)
    out = POST / f"{stem}_ens.per_gene_estimates.tsv"
    with open(out, "w") as fh:
        fh.write("ensg\tpost_mean\tprior_mean\tn_runs\n")
        for e in common:
            pm = sum(g[e] for g in gs) / len(gs)
            pr = sum(p[e] for p in ps) / len(ps)
            fh.write(f"{e}\t{pm:.10g}\t{pr:.10g}\t{len(runs)}\n")
    k = len(runs)
    rel = f"{sb(0.70, k):.2f}" if k > 1 else "0.70 (single run)"
    warn = "" if k >= 3 else "   <-- want >=3 runs for reliability >=0.87"
    print(f"[{stem}] averaged {k} converged run(s) over {len(common)} genes "
          f"-> {out.name}; est. reliability ~{rel}{warn}")
