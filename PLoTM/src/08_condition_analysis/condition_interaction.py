#!/usr/bin/env python
"""Rest-vs-Stim48hr condition interaction at the program level.

Programs are condition-specific (independent cNMF per condition), so we first MATCH programs across
conditions by reciprocal-best top-200 gene-set Jaccard. For each matched (rest_prog, stim_prog) pair
and each RA phenotype we compute the stimulation-interaction on the regulator-burden effect, which is
the only program-level quantity with a standard error:

    d_beta = beta_stim - beta_rest ,  SE_d = sqrt(SE_rest^2 + SE_stim^2) ,  z = d_beta/SE_d
    P_int  = 2*(1-Phi(|z|))                     BH-FDR within phenotype across matched pairs

This IS the G x condition contrast the request wants, expressed at the level the data actually vary
over (programs + beta), since gamma itself is condition-independent so a gene-level interaction is
identically zero and un-fittable.

Category per (pair, phenotype), using regulator FDR<0.1 as "present":
  constitutive        present both conditions, same sign, |d_beta| not FDR-sig
  stim-amplified      present both, same sign, d_beta FDR-sig and |beta_stim|>|beta_rest|
  stim-attenuated     present both, same sign, d_beta FDR-sig and |beta_stim|<|beta_rest|
  stim-revealed       present stim only
  rest-specific       present rest only
  sign-reversing      present both, opposite sign
  none                present in neither

Also meta-analyses d_beta across the 5 RA phenotypes (DL random effects) per matched pair.
"""
import numpy as np, pandas as pd
from scipy import stats

ROOT="/mnt/scratch/ZY2/Hackathon/pipeline_results"
CC=f"{ROOT}/condition_comparison"
RES="/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis"
tidy=pd.read_csv(f"{CC}/tables/tidy_results.tsv",sep="\t")
RA5=["BK_curated","BK_M06","BK_alt","GB_RA_custom","GB_RA_M06"]
SRC={"BK_curated":"Backman","BK_M06":"Backman","BK_alt":"Backman","GB_RA_custom":"GeneBass","GB_RA_M06":"GeneBass"}
TOP=200

# ---- reciprocal-best top-200 Jaccard match rest<->stim ----
def loadtop(ct):
    g=pd.read_csv(f"{RES}/cNMF/{ct}/test1/test1.gene_spectra_score.k_60.dt_0_4.txt",sep="\t",index_col=0)
    g.index=[f"P{i}" for i in g.index]
    return {p:set(g.loc[p].sort_values(ascending=False).index[:TOP]) for p in g.index}
R=loadtop("GWCD4i_rest_pseudobulk"); S=loadtop("GWCD4i_stim48_pseudobulk")
def jac(a,b): return len(a&b)/len(a|b)
J=pd.DataFrame({rp:{sp:jac(R[rp],S[sp]) for sp in S} for rp in R})  # rows=stim, cols=rest
pairs=[]
for rp in R:
    sp=J[rp].idxmax()                       # best stim for this rest
    if J.loc[sp].idxmax()==rp:              # reciprocal
        pairs.append((rp,sp,J.loc[sp,rp]))
pairs=[p for p in pairs if p[2]>=0.10]      # require modest overlap
lab=pd.read_csv("/mnt/scratch/ZY2/Hackathon/RA_cNMF_LoF_pipeline/program_annotation/curated_labels.tsv",sep="\t")
def ann(ct,p):
    m=lab[(lab.CT==ct)&(lab.Program==p)]; return m.Label.iloc[0] if len(m) else ""

def bh(p):
    p=np.asarray(p,float); m=np.isfinite(p); q=np.full(len(p),np.nan)
    pv=p[m]; n=len(pv); o=np.argsort(pv); r=pv[o]*n/np.arange(1,n+1)
    r=np.minimum.accumulate(r[::-1])[::-1]; qq=np.empty(n); qq[o]=np.clip(r,0,1); q[m]=qq; return q

def classify(brest,bstim,frest,fstim,dfdr):
    pr=np.isfinite(frest) and frest<0.1; ps=np.isfinite(fstim) and fstim<0.1
    if not pr and not ps: return "none"
    if pr and not ps: return "rest-specific"
    if ps and not pr: return "stim-revealed"
    if np.sign(brest)!=np.sign(bstim): return "sign-reversing"
    if np.isfinite(dfdr) and dfdr<0.1:
        return "stim-amplified" if abs(bstim)>abs(brest) else "stim-attenuated"
    return "constitutive"

rows=[]
for rp,sp,jj in pairs:
    for t in RA5:
        rr=tidy[(tidy.condition=="rest")&(tidy.short==t)&(tidy.program==rp)]
        ss=tidy[(tidy.condition=="stim48")&(tidy.short==t)&(tidy.program==sp)]
        if not len(rr) or not len(ss): continue
        rr=rr.iloc[0]; ss=ss.iloc[0]
        br,se_r=rr.reg_beta_withShet,rr.reg_betaSE_withShet
        bs,se_s=ss.reg_beta_withShet,ss.reg_betaSE_withShet
        se_d=np.sqrt(se_r**2+se_s**2); dbeta=bs-br; z=dbeta/se_d if se_d>0 else np.nan
        P=2*(1-stats.norm.cdf(abs(z))) if np.isfinite(z) else np.nan
        rows.append(dict(rest_prog=rp,stim_prog=sp,jaccard=jj,trait=t,source=SRC[t],
                         ann_rest=ann("GWCD4i_rest_pseudobulk",rp),ann_stim=ann("GWCD4i_stim48_pseudobulk",sp),
                         beta_rest=br,SE_rest=se_r,FDR_rest=rr.reg_FDR,
                         beta_stim=bs,SE_stim=se_s,FDR_stim=ss.reg_FDR,
                         d_beta=dbeta,SE_d=se_d,z=z,P_int=P))
D=pd.DataFrame(rows)
for t in RA5:
    m=D.trait==t; D.loc[m,"P_int_FDR"]=bh(D.loc[m,"P_int"].values)
D["category"]=[classify(r.beta_rest,r.beta_stim,r.FDR_rest,r.FDR_stim,r.P_int_FDR) for _,r in D.iterrows()]
D.to_csv(f"{CC}/tables/condition_interaction.tsv",sep="\t",index=False)

# meta of d_beta across the 5 traits per matched pair
def dl(beta,se):
    beta=np.asarray(beta,float);se=np.asarray(se,float);ok=np.isfinite(beta)&np.isfinite(se)&(se>0)
    beta,se=beta[ok],se[ok];k=len(beta)
    if k==0:return None
    w=1/se**2;fixed=np.sum(w*beta)/np.sum(w);Q=np.sum(w*(beta-fixed)**2);dfq=k-1
    C=np.sum(w)-np.sum(w**2)/np.sum(w) if k>1 else 0
    tau2=max(0,(Q-dfq)/C) if (k>1 and C>0) else 0.0
    I2=max(0,(Q-dfq)/Q)*100 if (k>1 and Q>0) else 0.0
    wr=1/(se**2+tau2);pooled=np.sum(wr*beta)/np.sum(wr);pse=np.sqrt(1/np.sum(wr))
    z=pooled/pse;p=2*(1-stats.norm.cdf(abs(z)))
    return dict(k=k,pooled=pooled,se=pse,P=p,I2=I2,Qp=(1-stats.chi2.cdf(Q,dfq)) if k>1 else np.nan)
mrows=[]
for (rp,sp),g in D.groupby(["rest_prog","stim_prog"]):
    m=dl(g.d_beta.values,g.SE_d.values)
    if m: mrows.append(dict(rest_prog=rp,stim_prog=sp,ann=g.ann_rest.iloc[0] or g.ann_stim.iloc[0],
                            jaccard=g.jaccard.iloc[0],pooled_dbeta=m["pooled"],SE=m["se"],P=m["P"],I2=m["I2"],Qp=m["Qp"],k=m["k"]))
MD=pd.DataFrame(mrows).sort_values("P")
MD["FDR"]=bh(MD.P.values)
MD.to_csv(f"{CC}/tables/condition_interaction_meta.tsv",sep="\t",index=False)

print(f"matched program pairs (reciprocal-best Jaccard>=0.10): {len(pairs)}")
print(f"condition_interaction.tsv: {len(D)} pair x trait rows")
print("\ncategory counts (per pair x trait):")
print(D.category.value_counts().to_string())
print("\nnon-'none' categories:")
print(D[D.category!="none"][["rest_prog","stim_prog","ann_rest","trait","source","beta_rest","beta_stim","d_beta","P_int_FDR","category"]].to_string(index=False))
print("\npooled d_beta (interaction) across 5 RA traits, top by P:")
print(MD.head(10)[["rest_prog","stim_prog","ann","pooled_dbeta","P","FDR","I2","k"]].to_string(index=False))
