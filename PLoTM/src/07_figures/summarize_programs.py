#!/usr/bin/env python
"""summarize_ra_programs.py — rank the cNMF programs by RA (rheumatoid-factor) burden and
annotate each top program with its top-loaded genes and its strongest regulators.

Consumes the verbatim Stage C outputs (program- and regulator-burden) plus the Stage A
gene_spectra_score and Stage B beta tables. Reporting only — computes no new statistics.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--ct", default="GWCD4i_stim48_pseudobulk")
ap.add_argument("--file", default="Genebass_RF.per_gene_estimates.tsv")
ap.add_argument("--topn-genes", type=int, default=12)
ap.add_argument("--topn-reg", type=int, default=10)
ap.add_argument("--show", type=int, default=8, help="how many top programs to detail")
args = ap.parse_args()

CT, FILE = args.ct, args.file
base = f"data/Perturbseq/trait_association/{CT}/ProgramLevel"
prog = pd.read_csv(f"{base}/programs_enrichment_K60_{FILE}", sep="\t")
reg = pd.read_csv(f"{base}/regulators_enrichment_K60_{FILE}", sep="\t")

# gencode ENSG<->symbol
corr = pd.read_csv("data/gencode_v41_gname_gid_ALL_sorted_onlyID", sep="\t", header=None,
                   names=["ensg", "sym"])
ensg2sym = dict(zip(corr.ensg, corr.sym))

# RF gamma
lof = pd.read_csv(f"data/LoF/GeneBayes_posterior/{FILE}", sep="\t")[["ensg", "post_mean"]]
ensg2gamma = dict(zip(lof.ensg, lof.post_mean))

# gene_spectra_score: rows=program 1..60, cols=ENSG genes
gs = pd.read_csv(
    f"data/Perturbseq/cNMF/{CT}/test1/test1.gene_spectra_score.k_60.dt_0_4.txt",
    sep="\t", index_col=0)
gs.index = [int(x) for x in gs.index]

# merge the two burden tables on Program
m = prog.merge(reg, on=["FILE", "Program"], how="outer").sort_values(
    "MEANgamma_top100_shet_adjusted_P")

def sym(e):
    return ensg2sym.get(e, e)

print("=" * 96)
print(f"RA (rheumatoid-factor) burden across {gs.shape[0]} cNMF programs — CT={CT}")
print(f"  program-burden  = mean RF gamma of a program's top-100 genes, shet-permuted P (100k)")
print(f"  regulator-burden= corr(beta_x(P), RF gamma) with shet-adjusted lm, P_withShet")
print("=" * 96)

cols = ["Program", "MEANgamma_top100", "shet_adjusted_random_mean",
        "MEANgamma_top100_shet_adjusted_P", "R_pearson_all", "P_withShet", "beta_withShet"]
disp = m[cols].copy()
disp.columns = ["Prog", "meanG_top100", "rand_mean", "progBurden_P", "reg_r", "reg_P_shet", "reg_beta"]
pd.set_option("display.width", 200, "display.max_columns", 20)
print("\n### All 60 programs, ranked by program-burden P (most RA-enriched first):\n")
print(disp.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

# Bonferroni threshold across 60 programs
bonf = 0.05 / gs.shape[0]
sig_prog = disp[disp.progBurden_P < bonf]
sig_reg = disp[disp.reg_P_shet < bonf]
print(f"\nBonferroni 0.05/60 = {bonf:.2e}")
print(f"  programs significant by PROGRAM-burden : {len(sig_prog)}  -> {list(sig_prog.Prog)}")
print(f"  programs significant by REGULATOR-burden: {len(sig_reg)}  -> {list(sig_reg.Prog)}")

# beta tables (Stage B)
def load_beta(p):
    return pd.read_csv(
        f"data/Perturbseq/cNMF_regulation/{CT}/K60_program{p}_perturb_effects.txt", sep="\t")

print("\n" + "=" * 96)
print(f"Top {args.show} programs by program-burden — top genes & strongest regulators")
print("=" * 96)
for _, row in m.head(args.show).iterrows():
    p = int(row["Program"])
    gp = float(row["MEANgamma_top100_shet_adjusted_P"])
    rp = float(row["P_withShet"]); rr = float(row["R_pearson_all"])
    # top-loaded genes
    srt = gs.loc[p].sort_values(ascending=False)
    top_g = [(sym(e), ensg2gamma.get(e, np.nan)) for e in srt.index[:args.topn_genes]]
    # strongest regulators (|beta|), annotate with RF gamma of the KD gene
    b = load_beta(p)
    b["absb"] = b.lm_es.abs()
    b = b.sort_values("absb", ascending=False).head(args.topn_reg)
    print(f"\n--- Program P{p} | progBurden_P={gp:.2e}  regBurden r={rr:+.3f} P_shet={rp:.2e} ---")
    print("  top-loaded genes (sym; RF gamma): " +
          ", ".join(f"{s}({g:+.3f})" if pd.notna(g) else f"{s}(NA)" for s, g in top_g))
    sym2ensg = {v: k for k, v in ensg2sym.items()}
    print("  strongest regulators (KD -> |beta|, lm_p, RF gamma of KD gene):")
    for _, br in b.iterrows():
        g = ensg2gamma.get(sym2ensg.get(br.GENE, ""), np.nan)
        gtxt = f"{g:+.3f}" if pd.notna(g) else "NA"
        print(f"      {br.GENE:12s} beta={br.lm_es:+.3f}  p={br.lm_p:.1e}  gamma={gtxt}")
print("\n[done]")
