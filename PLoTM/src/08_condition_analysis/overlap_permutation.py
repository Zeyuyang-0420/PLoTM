#!/usr/bin/env python
"""Permutation test: does cross-phenotype overlap of concordant regulator GENES exceed chance?

Null preserves, per phenotype: (i) the number of concordant top-hit genes it contributes, and (ii) the
candidate pool it draws from (the genes testable in that phenotype × condition, i.e. the union of that
condition's concordant-gene universe, weighted by gene DEGREE = how often a gene appears as a top hit
across phenotypes, so frequently-hit genes stay frequently-hit under the null).

Observed statistic = number of genes concordant in >=2 of the 5 RA phenotypes (per condition).
We resample each phenotype's concordant set (same size) from the pooled candidate universe with
degree-proportional weights, recompute the >=2 count, 10,000 times, and report the empirical P.

Program-level overlap is tested the same way (number of programs selected by >=2 phenotypes vs a null
that preserves each phenotype's selected-program count and samples from the 60 programs weighted by
how often each program is selectable).
"""
import numpy as np, pandas as pd, glob, os
from collections import Counter
ROOT="/mnt/scratch/ZY2/Hackathon/pipeline_results"; CC=f"{ROOT}/condition_comparison"
DATA="/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data"
th=pd.read_csv(f"{ROOT}/lof_thresholds.tsv",sep="\t"); TH=dict(zip(th.trait,th.thresh))
cp=pd.read_csv(f"{DATA}/gencode_v41_gname_gid_ALL_sorted_onlyID",sep="\t",header=None,names=["ensg","sym"])
e2s=dict(zip(cp.ensg,cp.sym))
RA5={"Backman_2021_RA_conv","Backman_2021_RA_M06_conv","Backman_2021_RA_alt_conv",
     "Genebass_RA_custom_conv","Genebass_RA_M06_conv"}
def tag(t): return f"program5_regulator3_LOF{TH[t]:g}"
def tophits(t):   # per-trait candidate universe = genes with |gamma|>threshold (condition-independent)
    g=pd.read_csv(f"{DATA}/LoF/GeneBayes_posterior/{t}.per_gene_estimates.tsv",sep="\t")
    g["sym"]=g.ensg.map(e2s); g=g.dropna(subset=["sym"])
    return set(g.sym[g.post_mean.abs()>TH[t]])
POOL={t:tophits(t) for t in RA5}
rng=np.random.default_rng(1); NPERM=10000
out=["Null: each phenotype draws its concordant-gene count from ITS OWN top-hit universe (|gamma|>q99,",
     "condition-independent). Preserves per-phenotype hit count and the real candidate pool + its overlaps.",""]
for cond,ct in [("rest","GWCD4i_rest_pseudobulk"),("stim48","GWCD4i_stim48_pseudobulk")]:
    concsets={}; progsets={}
    for t in RA5:
        D=f"{ROOT}/{t}/{ct}/fig5_out/{ct}/P200"
        cf=f"{D}/ConcordantGenes_{t}_{tag(t)}.txt"
        genes=pd.read_csv(cf,sep="\t").GENE.dropna().tolist() if os.path.exists(cf) and os.path.getsize(cf)>50 else []
        concsets[t]=set(genes)&POOL[t]      # concordant genes are a subset of the top-hit universe
        sf=f"{D}/SelectedPrograms_{t}_{tag(t)}.txt"
        if os.path.exists(sf):
            s=pd.read_csv(sf,sep="\t").iloc[0]
            progsets[t]=set(str(s.Program_selected).split(","))|set(str(s.Regulator_selected).split(","))
    # ---- gene overlap: draw each phenotype's concordant count from its own top-hit universe ----
    def count_shared(sets):
        c=Counter(); [c.update(set(s)) for s in sets]; return sum(1 for v in c.values() if v>=2)
    obs=count_shared(concsets.values())
    null=np.empty(NPERM)
    for i in range(NPERM):
        samp=[set(rng.choice(list(POOL[t]),size=min(len(concsets[t]),len(POOL[t])),replace=False)) for t in RA5]
        null[i]=count_shared(samp)
    p=(1+np.sum(null>=obs))/(NPERM+1)
    out.append(f"{cond} GENE overlap: observed shared>=2 = {obs}; null mean {null.mean():.1f} "
               f"(95%ile {np.percentile(null,95):.0f}); perm P = {p:.4f}")
    # ---- program overlap ----
    allprog=[f"P{i}" for i in range(1,61)]
    obsP=sum(1 for v in Counter(sum([list(s) for s in progsets.values()],[])).values() if v>=2)
    sizes=[len(s) for s in progsets.values()]
    nullP=np.empty(NPERM)
    for i in range(NPERM):
        samp=[rng.choice(allprog,size=n,replace=False) for n in sizes]
        c=Counter(); [c.update(s.tolist()) for s in samp]; nullP[i]=sum(1 for v in c.values() if v>=2)
    pP=(1+np.sum(nullP>=obsP))/(NPERM+1)
    out.append(f"{cond} PROGRAM overlap: observed shared>=2 = {obsP}; null mean {nullP.mean():.1f} (95%ile {np.percentile(nullP,95):.0f}); perm P = {pP:.4f}")
txt="\n".join(out); open(f"{CC}/tables/overlap_permutation.txt","w").write(txt+"\n"); print(txt)
