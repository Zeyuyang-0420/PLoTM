#!/usr/bin/env python
"""Harmonize the FULL unthresholded program-level burden results into one tidy table.

One row per (program, trait, condition). Columns:
  program, condition, trait, short, source, mask, program_class(RA-binary | RF-continuous), flag,
  annotation,
  # program-burden (membership route): mean gamma of top-100 program genes vs shet-matched null
  prog_meanG_top100, prog_null_mean, prog_effect(=meanG-null), prog_dir, prog_P, prog_FDR,
  # regulator-burden route: gamma ~ program-regulatory-beta (+ shet)
  reg_R_pearson, reg_P_pearson, reg_beta_withShet, reg_betaSE_withShet, reg_P_withShet, reg_dir, reg_FDR

FDR = Benjamini-Hochberg within each (trait, condition) across all 60 programs, separately for the
program-burden P and the regulator P.

NOTE on fields the request asks for that DO NOT EXIST in this pipeline's data (left out, not faked):
  gamma SE per gene exists (lower/upper_95) but is a gene-level LoF quantity, not per program-trait.
  carrier count / variant count / case-control counts: collapsed away upstream; absent from every file.
  burden mask: exactly ONE per trait (Backman M1_001, GeneBass M1_pLoF) -> no mask sensitivity.
  gamma is condition-independent -> Delta gamma == 0; condition signal is in programs+beta only.
"""
import os, glob
import numpy as np
import pandas as pd

ROOT = "/mnt/scratch/ZY2/Hackathon/pipeline_results"
LABF = "/mnt/scratch/ZY2/Hackathon/RA_cNMF_LoF_pipeline/program_annotation/curated_labels.tsv"
K = 60

# trait registry: canonical -> (short, source, mask, class, flag)
TRAITS = {
 "Genebass_RF_conv":                ("GB_RF",        "GeneBass","M1_pLoF","RF-continuous","ok"),
 "Genebass_RA_custom_conv":         ("GB_RA_custom", "GeneBass","M1_pLoF","RA-binary","ok"),
 "Genebass_RA_M06_conv":            ("GB_RA_M06",    "GeneBass","M1_pLoF","RA-binary","ok"),
 "Backman_2021_RA_conv":            ("BK_curated",   "Backman", "M1_001", "RA-binary","ok"),
 "Backman_2021_RA_M06_conv":        ("BK_M06",       "Backman", "M1_001", "RA-binary","BORDERLINE"),
 "Backman_2021_RA_alt_conv":        ("BK_alt",       "Backman", "M1_001", "RA-binary","ok"),
 "Backman_2021_RA_M05_seropos_conv":("BK_M05seropos","Backman", "M1_001", "RA-binary","PRIOR_DOMINATED"),
 "Backman_2021_RA_M060_seroneg_conv":("BK_M060seroneg","Backman","M1_001","RA-binary","SIGN_DEGENERATE"),
}
CONDS = [("GWCD4i_rest_pseudobulk","rest"), ("GWCD4i_stim48_pseudobulk","stim48")]

lab = pd.read_csv(LABF, sep="\t")
def ann(ct, prog):
    m = lab[(lab.CT==ct) & (lab.Program==f"P{prog}")]
    return m.Label.iloc[0] if len(m) else ""

def bh(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    ranked = p[order]*n/np.arange(1, n+1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(q, 0, 1); return out

rows = []
for trait,(short,source,mask,cls,flag) in TRAITS.items():
    for ct,cond in CONDS:
        PL = f"{ROOT}/{trait}/{ct}/data/Perturbseq/trait_association/{ct}/ProgramLevel"
        pf = f"{PL}/programs_enrichment_K{K}_{trait}.per_gene_estimates.tsv"
        rf = f"{PL}/regulators_enrichment_K{K}_{trait}.per_gene_estimates.tsv"
        if not (os.path.exists(pf) and os.path.exists(rf)):
            print(f"MISSING {trait}/{cond}"); continue
        p = pd.read_csv(pf, sep="\t"); r = pd.read_csv(rf, sep="\t")
        d = p.merge(r, on="Program", suffixes=("_p","_r"))
        d["prog_FDR"] = bh(d.MEANgamma_top100_shet_adjusted_P.values)
        d["reg_FDR"]  = bh(d.P_withShet.values)
        for _,x in d.iterrows():
            prog = int(x.Program)
            eff = x.MEANgamma_top100 - x.shet_adjusted_random_mean
            rows.append(dict(
                program=f"P{prog}", condition=cond, trait=trait, short=short, source=source,
                mask=mask, program_class=cls, flag=flag, annotation=ann(ct,prog),
                prog_meanG_top100=x.MEANgamma_top100, prog_null_mean=x.shet_adjusted_random_mean,
                prog_effect=eff, prog_dir=int(np.sign(eff)),
                prog_P=x.MEANgamma_top100_shet_adjusted_P, prog_FDR=x.prog_FDR,
                reg_R_pearson=x.R_pearson_all, reg_P_pearson=x.P_pearson_all,
                reg_beta_withShet=x.beta_withShet, reg_betaSE_withShet=x.betaSE_withShet,
                reg_P_withShet=x.P_withShet, reg_dir=int(np.sign(x.beta_withShet)), reg_FDR=x.reg_FDR))

df = pd.DataFrame(rows)
df.to_csv(f"{ROOT}/condition_comparison/tables/tidy_results.tsv", sep="\t", index=False)
print(f"tidy_results.tsv: {len(df)} rows = {df.trait.nunique()} traits x {df.condition.nunique()} conds x {K} programs")
print(f"  program-burden FDR<0.1: {(df.prog_FDR<0.1).sum()} rows ; regulator FDR<0.1: {(df.reg_FDR<0.1).sum()} rows")
print(df.groupby(['source','condition']).apply(lambda g: pd.Series({
    'prog_FDR<0.1': int((g.prog_FDR<0.1).sum()), 'reg_FDR<0.1': int((g.reg_FDR<0.1).sum())})).to_string())
