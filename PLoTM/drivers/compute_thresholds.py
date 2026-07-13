#!/usr/bin/env python
"""Per-trait |gamma| threshold + posterior QC for the 8 converged GeneBayes LoF posteriors.

Why per-trait: Backman gamma is ~60x the Genebass scale (sd 0.17-0.51 vs 0.005-0.008), so a fixed
|gamma| cutoff cannot be shared. |gamma|>0.03 keeps ~1% of genes for Genebass but 59-100% for
Backman -- for M060_seroneg it keeps ALL genes, leaving Fig 5a with an empty background.

We therefore use the 99th percentile of |gamma| per trait. This is the generalisation of the earlier
choice, not a departure from it: for Genebass_RF_conv the 99th percentile IS 0.0307.

QC flags (posterior vs the GeneBayes prior it shrinks toward):
  PRIOR_DOMINATED   r(post_mean, prior_mean) > 0.95  -> trait data barely move the posterior
  SIGN_DEGENERATE   <2% or >98% of gamma negative    -> sign(gamma) ~constant, so Fig 5a's
                                                        sign-concordance collapses into a
                                                        program-burden test
Idempotent: safe to re-run as new posteriors land.
"""
import os
import os, glob
import numpy as np
import pandas as pd

# Config-drivable via env (run_disease.sh sets these); falls back to the RA reference paths.
POST = os.environ.get("POSTERIOR_DIR",
    "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/RA_LoF/data/LoF/GeneBayes_posterior")
OUT = os.environ.get("THRESH_OUT",
    "/mnt/scratch/ZY2/Hackathon/pipeline_results/lof_thresholds.tsv")

# TRAITS: env-provided (comma or space list), else every *.per_gene_estimates.tsv in POST, else the RA panel.
_env = os.environ.get("TRAITS", "").replace(",", " ").split()
if _env:
    TRAITS = _env
elif glob.glob(f"{POST}/*.per_gene_estimates.tsv"):
    TRAITS = sorted(os.path.basename(p).replace(".per_gene_estimates.tsv", "")
                    for p in glob.glob(f"{POST}/*.per_gene_estimates.tsv"))
else:
    TRAITS = ["Genebass_RF_conv", "Genebass_RA_custom_conv", "Genebass_RA_M06_conv",
              "Backman_2021_RA_conv", "Backman_2021_RA_alt_conv", "Backman_2021_RA_M06_conv",
              "Backman_2021_RA_M05_seropos_conv", "Backman_2021_RA_M060_seroneg_conv"]

rows = []
for t in TRAITS:
    p = f"{POST}/{t}.per_gene_estimates.tsv"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        rows.append(dict(trait=t, status="PENDING", n=np.nan, thresh=np.nan,
                         n_top=np.nan, pct_neg=np.nan, r_prior=np.nan, flag=""))
        continue
    d = pd.read_csv(p, sep="\t")
    g = d.post_mean
    q = float(g.abs().quantile(0.99))
    pct_neg = round(100 * float((g < 0).mean()), 1)
    r_prior = round(float(np.corrcoef(g, d.prior_mean)[0, 1]), 3)
    f = []
    if pct_neg < 2 or pct_neg > 98:
        f.append("SIGN_DEGENERATE")
    if r_prior > 0.95:
        f.append("PRIOR_DOMINATED")
    # Backman_2021_RA_M06_conv lands just under both cutoffs (r_prior 0.937, 9.3% negative) while the
    # other unflagged traits sit at r_prior 0.89-0.92 and 16-48% negative. Not disqualifying, but it is
    # the most prior-influenced of the usable traits -- say so rather than lump it with the clean ones.
    if not f and (r_prior > 0.93 or pct_neg < 12):
        f.append("BORDERLINE")
    rows.append(dict(trait=t, status="READY", n=len(d), thresh=round(q, 5),
                     n_top=int((g.abs() > q).sum()), pct_neg=pct_neg, r_prior=r_prior,
                     flag=";".join(f) or "ok"))

df = pd.DataFrame(rows)
df["source"] = ["Backman" if t.startswith("Backman") else "Genebass" for t in df.trait]
# column order matters: run_batch.sh reads $1=trait $2=status $4=thresh
df = df[["trait", "status", "n", "thresh", "n_top", "pct_neg", "r_prior", "source", "flag"]]
df.to_csv(OUT, sep="\t", index=False)
pd.set_option("display.width", 220)
print(df.to_string(index=False, na_rep="."))
print(f"\nREADY {int((df.status=='READY').sum())}/8  -> {OUT}")
