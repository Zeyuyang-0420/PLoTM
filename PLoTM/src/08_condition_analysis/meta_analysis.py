#!/usr/bin/env python
"""Random-effects meta-analysis of the regulator-burden effect across the 5 binary RA phenotypes.

Effect = beta_withShet (gamma regressed on the program's regulatory-beta, shet-adjusted), with its SE
= betaSE_withShet. These are the ONLY program-level effects in this pipeline that carry a standard
error, so they are the only quantity a proper inverse-variance meta-analysis can use. (Program-burden
uses a permutation P with no SE; it is summarised by vote-counting / directional concordance instead.)

Per (program, condition) we meta-analyse across the 5 RA traits with DerSimonian-Laird random effects:
  pooled beta, SE, P ; Cochran Q + its P ; I^2 ; tau^2 ; leave-one-phenotype-out pooled betas.
GeneBass and Backman are NOT treated as independent cohorts: we report ALL-5, Backman-only (n=3),
GeneBass-only (n=2) and a cross-source agreement flag separately.

Excludes the two QC-failed Backman traits (M05_seropos PRIOR_DOMINATED, M060_seroneg SIGN_DEGENERATE).
RF is continuous and analysed on its own elsewhere, not pooled with the binary RA phenotypes.
BH-FDR is applied across programs within each condition on the ALL-5 pooled P.
"""
import numpy as np, pandas as pd
from scipy import stats

ROOT = "/mnt/scratch/ZY2/Hackathon/pipeline_results/condition_comparison"
df = pd.read_csv(f"{ROOT}/tables/tidy_results.tsv", sep="\t")

RA5 = ["BK_curated","BK_M06","BK_alt","GB_RA_custom","GB_RA_M06"]
SRC = {"BK_curated":"Backman","BK_M06":"Backman","BK_alt":"Backman",
       "GB_RA_custom":"GeneBass","GB_RA_M06":"GeneBass"}

def dl_meta(beta, se):
    """DerSimonian-Laird random-effects. Returns dict."""
    beta = np.asarray(beta,float); se = np.asarray(se,float)
    ok = np.isfinite(beta)&np.isfinite(se)&(se>0)
    beta,se = beta[ok],se[ok]; k=len(beta)
    if k==0: return None
    w = 1/se**2
    fixed = np.sum(w*beta)/np.sum(w)
    Q = np.sum(w*(beta-fixed)**2)
    df_q = k-1
    if k>1:
        C = np.sum(w)-np.sum(w**2)/np.sum(w)
        tau2 = max(0,(Q-df_q)/C) if C>0 else 0.0
        I2 = max(0,(Q-df_q)/Q)*100 if Q>0 else 0.0
        Qp = 1-stats.chi2.cdf(Q,df_q)
    else:
        tau2=0.0; I2=0.0; Qp=np.nan
    wr = 1/(se**2+tau2)
    pooled = np.sum(wr*beta)/np.sum(wr)
    pooled_se = np.sqrt(1/np.sum(wr))
    z = pooled/pooled_se; p = 2*(1-stats.norm.cdf(abs(z)))
    return dict(k=k, pooled=pooled, pooled_se=pooled_se, z=z, P=p, Q=Q, Q_P=Qp, I2=I2, tau2=tau2,
                n_pos=int(np.sum(beta>0)), n_neg=int(np.sum(beta<0)))

def bh(p):
    p=np.asarray(p,float); m=np.isfinite(p); q=np.full(len(p),np.nan)
    pv=p[m]; n=len(pv); o=np.argsort(pv); r=pv[o]*n/np.arange(1,n+1)
    r=np.minimum.accumulate(r[::-1])[::-1]; qq=np.empty(n); qq[o]=np.clip(r,0,1); q[m]=qq; return q

rows=[]
for cond in ["rest","stim48"]:
    sub = df[(df.condition==cond)&(df.short.isin(RA5))]
    for prog in [f"P{i}" for i in range(1,61)]:
        pr = sub[sub.program==prog]
        rec = dict(program=prog, condition=cond)
        for name,traits in [("all5",RA5),("Backman",[t for t in RA5 if SRC[t]=="Backman"]),
                            ("GeneBass",[t for t in RA5 if SRC[t]=="GeneBass"])]:
            g = pr[pr.short.isin(traits)]
            m = dl_meta(g.reg_beta_withShet.values, g.reg_betaSE_withShet.values)
            if m:
                for key in ["pooled","pooled_se","P","Q","Q_P","I2","tau2","n_pos","n_neg","k"]:
                    rec[f"{name}_{key}"]=m[key]
        # leave-one-phenotype-out on all5
        loo=[]
        for drop in RA5:
            g = pr[pr.short.isin([t for t in RA5 if t!=drop])]
            m = dl_meta(g.reg_beta_withShet.values, g.reg_betaSE_withShet.values)
            loo.append(m["pooled"] if m else np.nan)
        rec["loo_min"]=np.nanmin(loo); rec["loo_max"]=np.nanmax(loo)
        rec["loo_sign_stable"] = bool(np.all(np.sign([x for x in loo if np.isfinite(x)])==
                                             np.sign(rec.get("all5_pooled",np.nan)))) if np.isfinite(rec.get("all5_pooled",np.nan)) else False
        # cross-source agreement: Backman & GeneBass pooled betas same sign
        if np.isfinite(rec.get("Backman_pooled",np.nan)) and np.isfinite(rec.get("GeneBass_pooled",np.nan)):
            rec["cross_source_concordant"]= bool(np.sign(rec["Backman_pooled"])==np.sign(rec["GeneBass_pooled"]))
        else: rec["cross_source_concordant"]=False
        rows.append(rec)

M = pd.DataFrame(rows)
for cond in ["rest","stim48"]:
    mask = M.condition==cond
    M.loc[mask,"all5_FDR"] = bh(M.loc[mask,"all5_P"].values)
M = M.merge(df[["program","condition","annotation"]].drop_duplicates(), on=["program","condition"], how="left")
M.to_csv(f"{ROOT}/tables/meta_regulator_RA5.tsv", sep="\t", index=False)

# headline: significant pooled + heterogeneity
sig = M[(M.all5_FDR<0.1)].sort_values("all5_P")
print(f"meta_regulator_RA5.tsv: {len(M)} program x condition rows")
print(f"\nPooled regulator effect FDR<0.1 (all-5 RA):")
cols=["program","condition","annotation","all5_pooled","all5_P","all5_FDR","all5_I2","Q_p:=all5_Q_P",
      "cross_source_concordant","loo_sign_stable"]
if len(sig):
    print(sig[["program","condition","annotation","all5_pooled","all5_P","all5_FDR","all5_I2",
               "all5_Q_P","cross_source_concordant","loo_sign_stable"]].to_string(index=False))
else:
    print("  none at FDR<0.1")
print(f"\nlow-heterogeneity concordant programs (I2<25, cross-source concordant, |pooled|>0):")
lo = M[(M.all5_I2<25)&(M.cross_source_concordant)&(M.all5_pooled.abs()>0)].sort_values("all5_P").head(12)
print(lo[["program","condition","annotation","all5_pooled","all5_P","all5_I2","Backman_pooled","GeneBass_pooled"]].to_string(index=False))
