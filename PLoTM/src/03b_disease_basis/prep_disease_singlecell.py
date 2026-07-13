#!/usr/bin/env python
"""prep_disease_singlecell.py — build the cNMF input h5ad from disease single-cell data.

Reads either a disease .h5ad (raw counts in .X or a named layer) or the counts bundle exported by
rds_to_counts.R, optionally intersects genes with the Perturb-seq-measured set, applies the same
QC filter as Stage A filter_cells.py, and writes:

    data/perturbseq/filtered_data/<CT>.h5ad      (cells x genes, raw counts)  -- run_cnmf.sh input

so the disease basis is discovered by the SAME Stage A cNMF as everything else.

Usage:
    prep_disease_singlecell.py --config <cfg.yaml>
    prep_disease_singlecell.py --h5ad X.h5ad --ct RA_CD4T [--counts-layer counts] [--gene-id symbol]
                               [--perturb-h5ad perturbseq.h5ad] [--no-intersect]
"""
import argparse, os, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, anndata as ad
from scipy.sparse import csr_matrix
from scipy.io import mmread


def perturb_measured_genes(perturb_h5ad, gene_id):
    """Gene identifiers measured in the Perturb-seq atlas (symbols or ENSG), for the intersection."""
    a = sc.read_h5ad(perturb_h5ad, backed="r")
    if gene_id == "symbol" and "gene_name" in a.var.columns:
        return set(a.var["gene_name"].astype(str))
    return set(a.var_names.astype(str))


def load_disease(args):
    """Return an AnnData (cells x genes, raw counts) with obs carrying percent.mt/mitopercent if present."""
    prefix = args.rds_prefix
    if prefix and os.path.exists(prefix + ".counts.mtx"):
        M = mmread(prefix + ".counts.mtx").tocsr()                 # genes x cells
        genes = [l.strip() for l in open(prefix + ".genes.txt")]
        bcs = [l.strip() for l in open(prefix + ".barcodes.txt")]
        meta = pd.read_csv(prefix + ".metadata.csv").set_index("cell_barcode").loc[bcs]
        adata = ad.AnnData(X=csr_matrix(M.T), obs=meta,
                           var=pd.DataFrame(index=genes))
    else:
        adata = sc.read_h5ad(args.h5ad)
        if args.counts_layer and args.counts_layer in adata.layers:
            adata.X = adata.layers[args.counts_layer]
    adata.var_names_make_unique()
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config"); ap.add_argument("--h5ad"); ap.add_argument("--rds-prefix")
    ap.add_argument("--ct"); ap.add_argument("--counts-layer", default="counts")
    ap.add_argument("--gene-id", default="symbol", choices=["symbol", "ensg"])
    ap.add_argument("--perturb-h5ad"); ap.add_argument("--no-intersect", action="store_true")
    ap.add_argument("--data-dir", default="data")
    args = ap.parse_args()

    if args.config:
        import yaml
        cfg = yaml.safe_load(open(args.config))
        ds = cfg.get("disease_singlecell") or {}
        if not ds:
            sys.exit("[prep-sc] config has no disease_singlecell block — nothing to do")
        args.ct = args.ct or ds.get("ct_name", "disease")
        args.gene_id = ds.get("gene_id", args.gene_id)
        args.counts_layer = ds.get("counts_layer", args.counts_layer)
        args.h5ad = args.h5ad or ds.get("h5ad")
        args.perturb_h5ad = args.perturb_h5ad or cfg.get("perturbseq", {}).get("pseudobulk_h5ad")
        if not ds.get("intersect_perturb_measured", True):
            args.no_intersect = True
        if ds.get("rds") and not args.rds_prefix:
            args.rds_prefix = os.path.join(args.data_dir, "disease_sc", args.ct)

    adata = load_disease(args)
    print(f"[prep-sc] {args.ct}: raw {adata.shape} (cells x genes), gene_id={args.gene_id}", flush=True)

    # intersect with Perturb-seq-measured genes so the disease basis lives in a projectable gene space
    if not args.no_intersect and args.perturb_h5ad and os.path.exists(args.perturb_h5ad):
        pm = perturb_measured_genes(args.perturb_h5ad, args.gene_id)
        keep = [g for g in adata.var_names if g in pm]              # preserve disease gene order
        print(f"[prep-sc] intersect Perturb-measured: {len(keep)} / {adata.n_vars} genes", flush=True)
        adata = adata[:, keep].copy()

    # QC == Stage A filter_cells.py (mitopercent>0.3 removal, min_genes/counts=100, min_cells=10)
    if "mitopercent" not in adata.obs and "percent.mt" in adata.obs:
        adata.obs["mitopercent"] = adata.obs["percent.mt"].astype(float) / 100.0   # percent -> fraction
    if "mitopercent" in adata.obs:
        adata = adata[adata.obs["mitopercent"] <= 0.3].copy()
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_cells(adata, min_counts=100)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.X = csr_matrix(adata.X)
    print(f"[prep-sc] after QC: {adata.shape}; X integer={np.allclose(adata.X.data, np.round(adata.X.data))}", flush=True)

    out_dir = os.path.join(args.data_dir, "perturbseq", "filtered_data")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{args.ct}.h5ad")
    adata.write(out)
    print(f"[prep-sc] wrote {out}  -> now run: run_cnmf.sh {args.ct} <K>", flush=True)


if __name__ == "__main__":
    main()
