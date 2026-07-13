#!/usr/bin/env python
"""make_pseudobulk_input.py — build the cNMF/pipeline input from the GWCD4i pseudobulk h5ad.

Analogous to make_subset.py, but for the CD4+ T-cell PSEUDOBULK count matrix
(tcell_perturbseq/aggregate/GWCD4i.pseudobulk_merged.h5ad), which make_subset.py (single-cell)
cannot consume. We do NOT weaken any method: we only construct the input + the metadata schema the
verbatim pipeline scripts expect (cell_barcode, gene, gem_group, n_genes, mitopercent, leiden).

Scope (locked): keep_for_DE == True AND culture_condition == 'Stim48hr'.

Covariate mapping (pseudobulk -> the paper's single-cell schema):
  gene        = perturbed_gene_name  (NTC / guide_type 'non-targeting' -> 'non-targeting')
  gem_group   = donor_id + '_' + 10xrun_id           (batch; condition fixed at Stim48hr)
  n_genes     = # detected genes per profile (X>0)     (complexity/depth covariate)
  mitopercent = sum(MT- gene counts)/total_counts       (real QC covariate; 12 MT- genes in var)
  leiden      = 0  (carried, unused by any model)
"""
import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",
                    default="/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/aggregate/GWCD4i.pseudobulk_merged.h5ad")
    ap.add_argument("--ct", default="GWCD4i_stim48_pseudobulk")
    ap.add_argument("--condition", default="Stim48hr")
    args = ap.parse_args()

    print(f"[prep] reading {args.input} (backed) ...", flush=True)
    adata = sc.read_h5ad(args.input, backed="r")
    obs = adata.obs
    print(f"[prep] full: {adata.n_obs} profiles x {adata.n_vars} genes", flush=True)

    keep = (obs["keep_for_DE"].astype(bool).values) & (obs["culture_condition"].astype(str).values == args.condition)
    keep_idx = np.where(keep)[0]
    print(f"[prep] keep_for_DE & {args.condition}: {len(keep_idx)} profiles", flush=True)

    print("[prep] materializing subset into memory ...", flush=True)
    sub = adata[keep_idx].to_memory()
    if not sp.issparse(sub.X):
        sub.X = sp.csr_matrix(sub.X)
    sub.X = sub.X.tocsr()

    # drop the non-ENSG construct spike-in (CUSTOM001_PuroR); keep ENSG genes
    ensg_mask = np.array([str(v).startswith("ENSG") for v in sub.var_names])
    dropped = list(sub.var_names[~ensg_mask])
    sub = sub[:, ensg_mask].copy()
    sub.X = sub.X.tocsr()
    print(f"[prep] dropped {len(dropped)} non-ENSG vars {dropped}; now {sub.n_vars} genes", flush=True)

    counts = sub.X
    total = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()

    # mitopercent from MT- symbol genes
    gname = sub.var["gene_name"].astype(str)
    mt_mask = gname.str.upper().str.startswith("MT-").values
    mt_counts = np.asarray(counts[:, np.where(mt_mask)[0]].sum(axis=1)).ravel() if mt_mask.sum() else np.zeros_like(total)
    mitopercent = np.divide(mt_counts, total, out=np.zeros_like(total, dtype=float), where=total > 0)
    print(f"[prep] MT- genes used: {int(mt_mask.sum())}; "
          f"median n_genes={np.median(n_genes):.0f}, median mitopercent={np.median(mitopercent):.4f}, "
          f"max mitopercent={mitopercent.max():.4f}", flush=True)

    o = sub.obs
    gene = o["perturbed_gene_name"].astype(str).values.copy()
    is_ntc = (o["guide_type"].astype(str).values == "non-targeting")
    gene[is_ntc] = "non-targeting"
    gem_group = (o["donor_id"].astype(str) + "_" + o["10xrun_id"].astype(str)).values

    sub.obs = pd.DataFrame(
        {
            "gene": gene,
            "gem_group": gem_group,
            "n_genes": n_genes.astype(int),
            "mitopercent": mitopercent,
            "leiden": 0,
        },
        index=sub.obs_names.astype(str),
    )
    sub.obs.index.name = "cell_barcode"

    os.makedirs("data/Perturbseq/metadata", exist_ok=True)
    out_h5ad = f"data/Perturbseq/{args.ct}.h5ad"
    print(f"[prep] writing {out_h5ad}", flush=True)
    sub.write(out_h5ad)

    meta = sub.obs.copy()
    meta.index.name = "cell_barcode"
    out_meta = f"data/Perturbseq/metadata/{args.ct}_metadata.csv"
    print(f"[prep] writing {out_meta}", flush=True)
    meta.to_csv(out_meta)

    n_ctrl = int((sub.obs["gene"] == "non-targeting").sum())
    print("\n[prep] DONE", flush=True)
    print(f"  profiles       : {sub.n_obs}")
    print(f"  genes(var)     : {sub.n_vars}")
    print(f"  NTC controls   : {n_ctrl}")
    print(f"  perturbations  : {sub.obs['gene'].nunique() - 1}")
    print(f"  gem_groups     : {sub.obs['gem_group'].nunique()}  -> {sorted(sub.obs['gem_group'].unique())}")
    # how many perturbations pass Stage B's '>10 profiles' gate
    vc = sub.obs.loc[sub.obs['gene'] != 'non-targeting', 'gene'].value_counts()
    print(f"  perts with >10 profiles: {int((vc > 10).sum())}")


if __name__ == "__main__":
    main()
