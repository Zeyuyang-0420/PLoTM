#!/usr/bin/env python
"""Targeted priority-program / gene checks.

P31 IL-6/IL-1     : per-trait sign of the regulator + program burden in stim48; is GB_RA_M06 reversed?
P26, P1           : do GeneBass-selected programs stay directionally consistent in Backman BELOW the
                    selection threshold (full unthresholded burden, all 5 RA traits)?
P57 histone/mixed : top gene loadings + correlation of P57 usage with the cell-cycle program (P4 stim)
                    -> is P57 a cell-cycle/technical axis?
P58 olfactory     : correlate P58 usage with QC covariates available (n_genes = library complexity,
                    mitopercent = viability proxy, gem_group = batch). Flag for exclusion if QC-driven.
IFFO2 + shared    : gamma, cross-trait concordance; note infeasible variant/carrier LOO.
"""
import numpy as np, pandas as pd
from scipy import stats

ROOT="/mnt/scratch/ZY2/Hackathon/pipeline_results"; CC=f"{ROOT}/condition_comparison"
RES="/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis"
DATA="/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data"
tidy=pd.read_csv(f"{CC}/tables/tidy_results.tsv",sep="\t")
cp=pd.read_csv(f"{DATA}/gencode_v41_gname_gid_ALL_sorted_onlyID",sep="\t",header=None,names=["ensg","sym"])
e2s=dict(zip(cp.ensg,cp.sym))
RA5=["BK_curated","BK_M06","BK_alt","GB_RA_custom","GB_RA_M06"]
out=[]

def sec(t): out.append("\n"+"="*90+"\n"+t+"\n"+"="*90)

# ---------- P31: sign per trait (stim48) ----------
sec("P31 IL-6/IL-1 (stim48) — per-trait regulator & program-burden sign; is GB_RA_M06 reversed?")
p31=tidy[(tidy.condition=="stim48")&(tidy.program=="P31")&(tidy.short.isin(RA5))]
out.append(p31[["short","source","reg_beta_withShet","reg_P_withShet","reg_dir",
                "prog_effect","prog_P","prog_dir"]].to_string(index=False))
signs=p31.set_index("short").reg_dir
out.append(f"\n  regulator-burden signs: {dict(signs)}")
out.append(f"  GB_RA_M06 reg sign = {signs.get('GB_RA_M06')}, others = {sorted(set(signs.drop('GB_RA_M06',errors='ignore')))}"
           f"  => {'REVERSED' if signs.get('GB_RA_M06')!=stats.mode([s for k,s in signs.items() if k!='GB_RA_M06'],keepdims=False).mode else 'concordant with majority'}")

# ---------- P26 cytotoxic CD4 (stim P26) & P1 ribosome (stim P1): Backman sub-threshold consistency ----------
sec("P26 cytotoxic CD4 & P1 ribosome (stim48) — full unthresholded burden, GeneBass vs Backman")
for prog,name in [("P26","Cytotoxic CD4"),("P1","Ribosome/translation")]:
    sub=tidy[(tidy.condition=="stim48")&(tidy.program==prog)&(tidy.short.isin(RA5))]
    out.append(f"\n{prog} {name}:")
    out.append(sub[["short","source","prog_effect","prog_P","prog_dir","reg_beta_withShet","reg_P_withShet","reg_dir"]].to_string(index=False))
    gb=sub[sub.source=="GeneBass"].prog_dir; bk=sub[sub.source=="Backman"].prog_dir
    out.append(f"  program-burden dir: GeneBass={list(gb)} Backman={list(bk)} "
               f"-> Backman same-sign as GeneBass majority: {int((bk==stats.mode(gb,keepdims=False).mode).sum())}/{len(bk)}")

# ---------- P57 histone/mixed: loadings + cell-cycle correlation ----------
sec("P57 histone/mixed (stim48) — top gene loadings + correlation with cell-cycle usage")
g=pd.read_csv(f"{RES}/cNMF/GWCD4i_stim48_pseudobulk/test1/test1.gene_spectra_score.k_60.dt_0_4.txt",sep="\t",index_col=0)
g.index=[f"P{i}" for i in g.index]
top57=[e2s.get(x,x) for x in g.loc["P57"].sort_values(ascending=False).index[:15]]
out.append(f"  P57 top15 genes: {', '.join(top57)}")
u=pd.read_csv(f"{RES}/cNMF/GWCD4i_stim48_pseudobulk/test1/test1.usages.k_60.dt_0_4.consensus.txt",sep="\t",index_col=0)
u.columns=[f"P{c}" for c in u.columns]; u=u.div(u.sum(1),axis=0)   # normalise per profile
# stim48 cell-cycle program: find the program whose top genes are cell-cycle (MKI67/TOP2A/CCNB) — else use annotation
cc_prog=None
for p in g.index:
    tg=set(e2s.get(x,x) for x in g.loc[p].sort_values(ascending=False).index[:30])
    if len({"MKI67","TOP2A","CCNB1","CCNB2","CDK1","UBE2C","CENPF","BIRC5"}&tg)>=3: cc_prog=p;break
out.append(f"  identified cell-cycle program (stim48): {cc_prog} "
           f"(top: {', '.join([e2s.get(x,x) for x in g.loc[cc_prog].sort_values(ascending=False).index[:6]])})" if cc_prog else "  no clear cell-cycle program found")
if cc_prog:
    r=stats.spearmanr(u["P57"],u[cc_prog])
    out.append(f"  Spearman(P57 usage, {cc_prog} cell-cycle usage) = {r.correlation:.3f} (P={r.pvalue:.2e})")

# ---------- P58 olfactory/noise: QC correlation ----------
sec("P58 olfactory/noise (stim48) — correlation of program usage with QC covariates")
meta=pd.read_csv(f"{RES}/GWCD4i_stim48_pseudobulk_metadata.csv").set_index("cell_barcode")
uj=u.join(meta[["n_genes","mitopercent","gem_group"]],how="inner")
for qc,lab in [("n_genes","library complexity (n_genes)"),("mitopercent","viability (mito%)")]:
    r=stats.spearmanr(uj["P58"],uj[qc]); out.append(f"  Spearman(P58 usage, {lab}) = {r.correlation:+.3f} (P={r.pvalue:.1e})")
# batch: eta^2 of P58 usage explained by gem_group (one-way ANOVA)
grps=[v["P58"].values for _,v in uj.groupby("gem_group")]
F,pF=stats.f_oneway(*grps)
ss_between=sum(len(gg)*(gg.mean()-uj["P58"].mean())**2 for gg in grps); ss_tot=((uj["P58"]-uj["P58"].mean())**2).sum()
out.append(f"  batch (gem_group) one-way ANOVA on P58 usage: F={F:.1f} P={pF:.1e}  eta^2={ss_between/ss_tot:.3f}")
# compare to a real biological program (P26 cytotoxic) as reference
for ref in ["P26"]:
    grp2=[v[ref].values for _,v in uj.groupby("gem_group")]
    ssb=sum(len(gg)*(gg.mean()-uj[ref].mean())**2 for gg in grp2); sst=((uj[ref]-uj[ref].mean())**2).sum()
    rr=stats.spearmanr(uj[ref],uj["n_genes"])
    out.append(f"  [ref {ref}] Spearman(usage,n_genes)={rr.correlation:+.3f}; batch eta^2={ssb/sst:.3f}")
top58=[e2s.get(x,x) for x in g.loc["P58"].sort_values(ascending=False).index[:12]]
out.append(f"  P58 top12 genes: {', '.join(top58)}")

# ---------- IFFO2 + shared regulator genes ----------
sec("IFFO2 + shared regulator genes — feasible checks (variant/carrier LOO NOT possible: no such data)")
for gene in ["IFFO2","ATXN1","NFIA"]:
    ens=[k for k,v in e2s.items() if v==gene]
    line=f"  {gene}: "
    for t in ["Genebass_RA_M06_conv","Genebass_RA_custom_conv","Backman_2021_RA_conv"]:
        try:
            gg=pd.read_csv(f"{DATA}/LoF/GeneBayes_posterior/{t}.per_gene_estimates.tsv",sep="\t")
            row=gg[gg.ensg.isin(ens)]
            if len(row): line+=f"{t.replace('_conv','')} gamma={row.post_mean.iloc[0]:+.4f}; "
        except Exception: pass
    out.append(line)
out.append("  NOTE: leave-one-variant-out / leave-one-carrier-out require per-variant / per-carrier data")
out.append("        that GeneBayes collapsed away; only gene-level gamma (+/-95% CI) survives. Not feasible here.")

txt="\n".join(out)
open(f"{CC}/tables/priority_checks.txt","w").write(txt)
print(txt)
