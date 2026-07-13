#!/usr/bin/env python
"""Figures for the Rest-vs-Stim48hr condition comparison.

1 program_heatmap        programs x (trait,condition); tile = standardized regulator beta; * = FDR<0.1;
                         source encoded on the trait axis; side bars = recurrence, concordance, I2.
2 rest_vs_stim_scatter   per RA trait: x = signed program-burden -log10P Rest, y = Stim48hr (matched
                         programs by Jaccard); identity line; label large |delta| / FDR interactions.
3 forest_priority        priority programs: per-trait regulator beta +/-95%CI + DL pooled (all5 / by source).
4 upset_sign             same-sign vs opposite-sign overlap of nominal regulator hits across the 5 traits.
5 p58_qc                 P58 usage vs mito%, n_genes, batch — is it QC-driven?
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy import stats

CC="/mnt/scratch/ZY2/Hackathon/pipeline_results/condition_comparison"
RES="/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis"
FIG=f"{CC}/figures"
tidy=pd.read_csv(f"{CC}/tables/tidy_results.tsv",sep="\t")
tidy["annotation"]=tidy["annotation"].fillna("")   # empty labels round-trip to NaN; treat as ""
meta=pd.read_csv(f"{CC}/tables/meta_regulator_RA5.tsv",sep="\t")
RA5=["BK_curated","BK_M06","BK_alt","GB_RA_custom","GB_RA_M06"]
SRC={"BK_curated":"Backman","BK_M06":"Backman","BK_alt":"Backman","GB_RA_custom":"GeneBass","GB_RA_M06":"GeneBass"}

# ---------------- 1. program heatmap ----------------
def heatmap():
    # TWO INDEPENDENT per-condition panels. cNMF ran separately per condition, so rest-Pi and stim-Pi are
    # DIFFERENT programs and must NOT share a row. Each panel is ranked & annotated WITHIN its condition;
    # cross-condition comparison is valid ONLY via the Jaccard-matched pairs in fig 1b.
    def panel_data(cond):
        sub=tidy[tidy.condition==cond]
        annp=set(sub[sub.annotation!=""].program)                     # annotated IN THIS condition
        keep=sorted(sub.program.unique(),key=lambda p:(p not in annp,int(p[1:])))  # ALL 60, annotated first
        Z=np.full((len(keep),len(RA5)),np.nan); FDR=np.full_like(Z,np.nan); ylab=[]
        for i,p in enumerate(keep):
            a=sub[(sub.program==p)&(sub.annotation!="")].annotation
            ylab.append(f"{p}  {a.iloc[0]}" if len(a) else p)
            for j,t in enumerate(RA5):
                r=sub[(sub.program==p)&(sub.short==t)]
                if len(r): Z[i,j]=r.reg_beta_withShet.iloc[0]; FDR[i,j]=r.reg_FDR.iloc[0]
        Zs=(Z-np.nanmean(Z,0))/np.nanstd(Z,0)                          # standardize per phenotype column
        return keep,ylab,Zs,FDR
    P={c:panel_data(c) for c in ["rest","stim48"]}
    nmax=max(len(P["rest"][0]),len(P["stim48"][0]))
    vmax=np.nanpercentile(np.abs(np.concatenate([P["rest"][2].ravel(),P["stim48"][2].ravel()])),98)
    fig,axes=plt.subplots(1,2,figsize=(15,max(6,nmax*0.30)),gridspec_kw=dict(wspace=0.6))
    for ax,cond,title in [(axes[0],"rest","REST"),(axes[1],"stim48","STIM48hr")]:
        keep,ylab,Zs,FDR=P[cond]
        im=ax.imshow(Zs,aspect="auto",cmap="RdBu_r",vmin=-vmax,vmax=vmax)
        ax.set_xticks(range(len(RA5))); ax.set_xticklabels(RA5,rotation=90,fontsize=8)
        ax.set_yticks(range(len(keep))); ax.set_yticklabels(ylab,fontsize=7)
        for i in range(len(keep)):
            for j in range(len(RA5)):
                if np.isfinite(FDR[i,j]) and FDR[i,j]<0.1:
                    ax.text(j,i,"*",ha="center",va="center",fontsize=10,fontweight="bold")
        for j,t in enumerate(RA5):
            c="#4c72b0" if SRC[t]=="GeneBass" else "#c44e52"
            ax.add_patch(Rectangle((j-0.5,-0.9),1,0.4,color=c,clip_on=False))
        ax.set_title(f"{title}\n(programs ranked within this condition)",fontsize=10,fontweight="bold")
    cb=fig.colorbar(im,ax=axes,fraction=0.02,pad=0.02); cb.set_label("standardized regulator β (per phenotype column)",fontsize=8)
    fig.suptitle("Regulator-burden per program × phenotype, EACH CONDITION SEPARATELY\n"
                 "rows independent per panel (cNMF run separately) — do NOT read across panels; use 1b for matched programs.\n"
                 "* = FDR<0.1 ; blue x-strip = GeneBass, red = Backman ; β>0 = regulator-LoF↑ → program↑ → trait↑",fontsize=9.5,y=1.03)
    fig.savefig(f"{FIG}/1_program_heatmap.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/1_program_heatmap.pdf",bbox_inches="tight"); plt.close(fig)
    print("  1_program_heatmap")

# ---------------- 2. rest vs stim scatter (per trait) ----------------
def scatter():
    ci=pd.read_csv(f"{CC}/tables/condition_interaction.tsv",sep="\t")
    fig,axes=plt.subplots(1,5,figsize=(22,4.6),sharex=True,sharey=True)
    for ax,t in zip(axes,RA5):
        g=ci[ci.trait==t]
        x=g.beta_rest; y=g.beta_stim
        lim=max(abs(pd.concat([x,y])).max()*1.1,0.02)
        ax.axhline(0,color="grey",lw=.5); ax.axvline(0,color="grey",lw=.5)
        ax.plot([-lim,lim],[-lim,lim],"--",color="grey",lw=.8)
        sig=g.P_int_FDR<0.1
        ax.scatter(x[~sig],y[~sig],s=22,c="#999",alpha=.6,edgecolor="none")
        ax.scatter(x[sig],y[sig],s=55,c="#d1495b",edgecolor="k",zorder=3)
        for _,r in g[sig | (g.d_beta.abs()>g.d_beta.abs().quantile(.9))].iterrows():
            lab=r.ann_rest if isinstance(r.ann_rest,str) and r.ann_rest else r.rest_prog
            ax.annotate(f"{r.rest_prog}/{r.stim_prog} {str(lab)[:14]}",(r.beta_rest,r.beta_stim),
                        fontsize=6,xytext=(3,3),textcoords="offset points")
        ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
        ax.set_title(f"{t} ({SRC[t]})",fontsize=10); ax.set_xlabel("regulator β  Rest")
    axes[0].set_ylabel("regulator β  Stim48hr")
    fig.suptitle("Rest vs Stim48hr regulator-burden β per matched program (red = FDR<0.1 interaction; dashed = identity)",fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/2_rest_vs_stim_scatter.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/2_rest_vs_stim_scatter.pdf",bbox_inches="tight"); plt.close(fig)
    print("  2_rest_vs_stim_scatter")

# ---------------- 3. forest for priority programs ----------------
def forest():
    prio=[("rest","P2","Translation init"),("rest","P24","TNF-NFkB"),("rest","P45","Myeloid S100A9"),
          ("rest","P35","IFN-alpha"),("rest","P49","Prostaglandin-RANKL-Th2"),
          ("stim48","P31","IL-6/IL-1"),("stim48","P26","Cytotoxic CD4"),("stim48","P1","Ribosome")]
    fig,axes=plt.subplots(1,len(prio),figsize=(24,4.4),sharey=False)
    for ax,(cond,prog,name) in zip(axes,prio):
        rows=[]
        for t in RA5:
            r=tidy[(tidy.condition==cond)&(tidy.program==prog)&(tidy.short==t)]
            if len(r): rows.append((t,SRC[t],r.reg_beta_withShet.iloc[0],r.reg_betaSE_withShet.iloc[0]))
        rows=rows[::-1]
        y=np.arange(len(rows))
        for i,(t,s,b,se) in enumerate(rows):
            c="#4c72b0" if s=="GeneBass" else "#c44e52"
            ax.errorbar(b,i,xerr=1.96*se,fmt="o",color=c,capsize=2,ms=5)
            ax.text(b,i+0.18,t,fontsize=6,ha="center")
        m=meta[(meta.condition==cond)&(meta.program==prog)]
        if len(m):
            mp=m.iloc[0]
            for name2,key,yy,col in [("all5","all5",-1.3,"k"),("Backman","Backman",-2.0,"#c44e52"),("GeneBass","GeneBass",-2.7,"#4c72b0")]:
                b=mp.get(f"{key}_pooled",np.nan); se=mp.get(f"{key}_pooled_se",np.nan)
                if np.isfinite(b) and np.isfinite(se):
                    ax.errorbar(b,yy,xerr=1.96*se,fmt="D",color=col,capsize=3,ms=6)
                    ax.text(b,yy-0.35,f"{name2}",fontsize=6,ha="center",color=col)
        ax.axvline(0,color="grey",lw=.6)
        i2=meta[(meta.condition==cond)&(meta.program==prog)].all5_I2
        ax.set_title(f"{prog} {name}\n({cond}) I²={i2.iloc[0]:.0f}%" if len(i2) else f"{prog} {name}",fontsize=8)
        ax.set_yticks([]); ax.set_xlabel("regulator β",fontsize=8); ax.tick_params(labelsize=7)
    fig.suptitle("Priority programs: per-phenotype regulator β ±95% CI (circles) and DL pooled (diamonds)",fontsize=11)
    fig.tight_layout(); fig.savefig(f"{FIG}/3_forest_priority.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/3_forest_priority.pdf",bbox_inches="tight"); plt.close(fig)
    print("  3_forest_priority")

# ---------------- 4. sign-overlap (same vs opposite) across traits ----------------
def upset_sign():
    fig,axes=plt.subplots(1,2,figsize=(13,5))
    for ax,cond in zip(axes,["rest","stim48"]):
        sub=tidy[(tidy.condition==cond)&(tidy.short.isin(RA5))]
        progs=[f"P{i}" for i in range(1,61)]
        # for each program: n traits with nominal reg P<0.05, and sign agreement
        rec=[]
        for p in progs:
            g=sub[sub.program==p]
            hit=g[g.reg_P_withShet<0.05]
            if len(hit)>=2:
                npos=(hit.reg_dir>0).sum(); nneg=(hit.reg_dir<0).sum()
                rec.append((p,g.annotation.iloc[0],len(hit),npos,nneg))
        rec=sorted(rec,key=lambda x:-x[2])[:14]
        y=np.arange(len(rec))
        ax.barh(y,[r[3] for r in rec],color="#c0392b",label="β>0 (LoF↑→trait↑)")
        ax.barh(y,[-r[4] for r in rec],color="#2c6fbb",label="β<0")
        ax.set_yticks(y); ax.set_yticklabels([f"{r[0]} {(str(r[1])[:16] if isinstance(r[1],str) else '')}" for r in rec],fontsize=7)
        ax.axvline(0,color="k",lw=.6); ax.set_title(f"{cond}: programs with ≥2 nominal hits (n traits by sign)",fontsize=9)
        ax.set_xlabel("# RA phenotypes (nominal reg P<0.05)"); ax.legend(fontsize=7,loc="lower right")
    fig.tight_layout(); fig.savefig(f"{FIG}/4_sign_overlap.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/4_sign_overlap.pdf",bbox_inches="tight"); plt.close(fig)
    print("  4_sign_overlap")

# ---------------- 5. P58 QC ----------------
def p58qc():
    g=pd.read_csv(f"{RES}/cNMF/GWCD4i_stim48_pseudobulk/test1/test1.gene_spectra_score.k_60.dt_0_4.txt",sep="\t",index_col=0)
    u=pd.read_csv(f"{RES}/cNMF/GWCD4i_stim48_pseudobulk/test1/test1.usages.k_60.dt_0_4.consensus.txt",sep="\t",index_col=0)
    u.columns=[f"P{c}" for c in u.columns]; u=u.div(u.sum(1),axis=0)
    md=pd.read_csv(f"{RES}/GWCD4i_stim48_pseudobulk_metadata.csv").set_index("cell_barcode")
    d=u.join(md[["n_genes","mitopercent","gem_group"]],how="inner")
    fig,axes=plt.subplots(1,3,figsize=(15,4.4))
    for ax,(qc,lab) in zip(axes,[("mitopercent","mito fraction (viability)"),("n_genes","n_genes (complexity)")]):
        ax.scatter(d[qc],d["P58"],s=3,alpha=.15,c="#666")
        r=stats.spearmanr(d[qc],d["P58"])
        ax.set_xlabel(lab); ax.set_ylabel("P58 usage (normalised)")
        ax.set_title(f"P58 vs {lab}\nSpearman ρ={r.correlation:+.3f}",fontsize=9)
    # batch boxplot
    ax=axes[2]
    order=sorted(d.gem_group.unique())
    ax.boxplot([d[d.gem_group==b]["P58"].values for b in order],showfliers=False)
    ax.set_xticklabels(order,rotation=90,fontsize=6); ax.set_ylabel("P58 usage"); ax.set_title("P58 usage by batch (gem_group)",fontsize=9)
    fig.suptitle("P58 (olfactory-receptor genes) QC dependence — top genes: OR10G2, OR8A1, ... → exclude from biology",fontsize=10)
    fig.tight_layout(); fig.savefig(f"{FIG}/5_p58_qc.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/5_p58_qc.pdf",bbox_inches="tight"); plt.close(fig)
    print("  5_p58_qc")

# ---------------- 1b. condition-CHANGE heatmap (Rest -> Stim48hr) ----------------
def change_heatmap():
    ci=pd.read_csv(f"{CC}/tables/condition_interaction.tsv",sep="\t")
    # one row per matched program pair; annotate; annotated pairs first
    pairs=ci[["rest_prog","stim_prog","jaccard","ann_rest","ann_stim"]].drop_duplicates()
    def plabel(r):
        a=r.ann_rest if isinstance(r.ann_rest,str) and r.ann_rest else (r.ann_stim if isinstance(r.ann_stim,str) else "")
        return f"{r.rest_prog}/{r.stim_prog}"+(f"  {a}" if a else "")
    pairs["label"]=pairs.apply(plabel,axis=1)
    pairs["annotated"]=pairs.apply(lambda r:bool((isinstance(r.ann_rest,str) and r.ann_rest) or (isinstance(r.ann_stim,str) and r.ann_stim)),axis=1)
    pairs=pairs.sort_values(["annotated","jaccard"],ascending=[False,False]).reset_index(drop=True)
    # matrix of d_beta (stim - rest) and interaction significance
    D=np.full((len(pairs),len(RA5)),np.nan); Pnom=np.full_like(D,np.nan); Pfdr=np.full_like(D,np.nan)
    for i,pr in pairs.iterrows():
        for j,t in enumerate(RA5):
            r=ci[(ci.rest_prog==pr.rest_prog)&(ci.stim_prog==pr.stim_prog)&(ci.trait==t)]
            if len(r): D[i,j]=r.d_beta.iloc[0]; Pnom[i,j]=r.P_int.iloc[0]; Pfdr[i,j]=r.P_int_FDR.iloc[0]
    fig,ax=plt.subplots(figsize=(9.5,max(6,len(pairs)*0.34)))
    vmax=np.nanpercentile(np.abs(D),98)
    im=ax.imshow(D,aspect="auto",cmap="PuOr_r",vmin=-vmax,vmax=vmax)   # orange=up in stim, purple=down
    # x labels carry the per-column change summary: #nominal-sig up / down
    xlab=[]
    for j,t in enumerate(RA5):
        col=D[:,j]; sig=(Pnom[:,j]<0.05)
        nUp=int(np.nansum((col>0)&sig)); nDn=int(np.nansum((col<0)&sig))
        xlab.append(f"{t}\n▲{nUp} ▼{nDn}")
    ax.set_xticks(range(len(RA5))); ax.set_xticklabels(xlab,rotation=0,fontsize=7.5)
    ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs.label,fontsize=7)
    # on each tile: direction arrow, doubled/starred by significance
    for i in range(len(pairs)):
        for j in range(len(RA5)):
            if not np.isfinite(D[i,j]): continue
            arr="▲" if D[i,j]>0 else "▼"                      # higher in Stim / higher in Rest
            if Pfdr[i,j]<0.1:  ax.text(j,i,f"{arr}{arr}",ha="center",va="center",fontsize=9,fontweight="bold")
            elif Pnom[i,j]<0.05: ax.text(j,i,f"{arr}*",ha="center",va="center",fontsize=8,fontweight="bold")
    # source strip just above the x labels
    for j,t in enumerate(RA5):
        c="#4c72b0" if SRC[t]=="GeneBass" else "#c44e52"
        ax.add_patch(Rectangle((j-0.5,len(pairs)-0.5),1,0.30,color=c,clip_on=False))
    cb=fig.colorbar(im,ax=ax,fraction=0.03,pad=0.02)
    cb.set_label("Δβ = β(Stim48hr) − β(Rest)   orange ▲ higher in Stim ·  purple ▼ higher in Rest",fontsize=8)
    ax.set_title("Rest → Stim48hr CHANGE in regulator-burden β, per matched program × phenotype\n"
                 "tile arrow = direction (▲ up in Stim, ▼ up in Rest); ▲▲ interaction FDR<0.1, ▲* nominal P<0.05\n"
                 "rows = Jaccard-matched pairs Prest/Pstim (annotated at top) · x-strip blue=GeneBass red=Backman · ▲/▼ counts under each column = # nominal-sig",
                 fontsize=8.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/1b_condition_change_heatmap.png",dpi=150,bbox_inches="tight")
    fig.savefig(f"{FIG}/1b_condition_change_heatmap.pdf",bbox_inches="tight"); plt.close(fig)
    print("  1b_condition_change_heatmap")

heatmap(); change_heatmap(); scatter(); forest(); upset_sign(); p58qc()
print("all figures ->",FIG)
