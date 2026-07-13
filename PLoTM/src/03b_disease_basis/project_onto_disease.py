"""project_onto_disease.py — project a Perturb-seq condition onto the FIXED disease cNMF basis.

The disease basis W = the disease cNMF's `gene_spectra_tpm` (K x genes, non-negative, TPM units).
For a Perturb-seq condition, X = TPM(condition counts)[:, disease genes] (library-size normalised,
no log, no scaling — exactly how cNMF built the matrix its usages were solved against), and

    U = argmin_{U>=0} || X - U W ||_F^2      (sklearn NNLS; W fixed, update_H=False)

W is NEVER updated (assert-checked). Then, so Stage B / Stage C / figures run UNCHANGED, we write
into the condition's CT path the same file set the native cNMF stage produces:

    <RESULTS>/cNMF/<CT>/test1/test1.usages.k_<K>.dt_0_4.consensus.txt      (profiles x K)  <- U
    <RESULTS>/cNMF/<CT>/test1/test1.gene_spectra_score.k_<K>.dt_0_4.txt    (K x genes)     <- disease's
    <RESULTS>/cNMF/<CT>/test1/test1.gene_spectra_tpm.k_<K>.dt_0_4.txt      (K x genes)     <- disease's
    <RESULTS>/cNMF/<CT>/program_estimability.tsv                          (per program NTC usage stats)

Program identity (gene_spectra_score) is the DISEASE's, shared across conditions; only the usages
differ per condition. Non-estimable programs (zero usage variance among NTC controls) have no
Stage-B beta and are flagged here so downstream can skip them.

Usage:
    project_onto_disease.py --disease-ct RA_CD4T --disease-k 30 \
        --condition-h5ad data/perturbseq/filtered_data/<CT>.h5ad --ct <CT> \
        --metadata data/Perturbseq/metadata/<CT>_metadata.csv [--results <RESULTS_ROOT>]
"""
import argparse, os, warnings, yaml
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from sklearn.decomposition import non_negative_factorization


def nnls_fixed_basis(X, W, nmf_kwargs, seed=14):
    kw = dict(nmf_kwargs)
    kw.update(dict(n_components=W.shape[0], H=np.ascontiguousarray(W, dtype=X.dtype),
                   update_H=False, random_state=seed))
    U, H, _ = non_negative_factorization(X, **kw)
    assert np.allclose(H, W, atol=1e-10), "disease basis W was modified — it must stay fixed"
    return U


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease-ct", required=True)
    ap.add_argument("--disease-k", type=int, required=True)
    ap.add_argument("--condition-h5ad", required=True)
    ap.add_argument("--ct", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--results", default=os.environ.get("RESULTS_ROOT", "results"))
    args = ap.parse_args()
    K = args.disease_k
    dis = os.path.join(args.results, "cNMF", args.disease_ct, "test1")

    # --- disease basis: gene_spectra_tpm (projection) + tpm_stats (unused here) ---
    W_full = pd.read_csv(os.path.join(dis, f"test1.gene_spectra_tpm.k_{K}.dt_0_4.txt"),
                         sep="\t", index_col=0)                     # K x genes, TPM units, >=0
    genes = list(W_full.columns)
    score = pd.read_csv(os.path.join(dis, f"test1.gene_spectra_score.k_{K}.dt_0_4.txt"),
                        sep="\t", index_col=0)                      # K x genes, program identity
    nmf_kwargs = yaml.safe_load(open(os.path.join(dis, "cnmf_tmp", "test1.nmf_idvrun_params.yaml")))
    print(f"[proj] disease basis {args.disease_ct}: W {W_full.shape} (K={K})", flush=True)

    # --- condition counts -> TPM on the shared disease gene universe ---
    a = sc.read_h5ad(args.condition_h5ad)
    if "gene_name" in a.var.columns:                                # ENSG-indexed Perturb h5ad -> symbols
        a.var_names = a.var["gene_name"].astype(str).values
    common = [g for g in genes if g in set(a.var_names)]            # keep W's gene order
    missing = [g for g in genes if g not in set(a.var_names)]
    print(f"[proj] {args.ct}: common with disease gene space {len(common)}/{len(genes)}"
          + (f"; {len(missing)} disease genes absent -> dropped from BOTH X and W" if missing else ""), flush=True)
    W = W_full.loc[:, common]
    a = a[:, common].copy()
    sc.pp.normalize_total(a, target_sum=1e6)                        # == cnmf.compute_tpm
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    X = X.astype(np.float64)
    cells = a.obs_names.to_list(); del a

    U = nnls_fixed_basis(X, W.values, nmf_kwargs)
    R = X - U @ W.values
    rel = np.sqrt(np.einsum("ij,ij->i", R, R) / np.maximum(np.einsum("ij,ij->i", X, X), 1e-300))
    dead = [int(W.index[k]) for k in range(K) if (U[:, k] > 0).sum() == 0]
    print(f"[proj] U {U.shape} all>=0 rel_err median={np.median(rel):.4f} dead_programs={dead}", flush=True)

    # --- write the CT file set the downstream expects (usages + program identity) ---
    out = os.path.join(args.results, "cNMF", args.ct, "test1")
    os.makedirs(os.path.join(out, "cnmf_tmp"), exist_ok=True)
    U_df = pd.DataFrame(U, index=cells, columns=[str(c) for c in W.index]); U_df.index.name = "cell_barcode"
    U_df.to_csv(os.path.join(out, f"test1.usages.k_{K}.dt_0_4.consensus.txt"), sep="\t")
    score.to_csv(os.path.join(out, f"test1.gene_spectra_score.k_{K}.dt_0_4.txt"), sep="\t")
    W_full.to_csv(os.path.join(out, f"test1.gene_spectra_tpm.k_{K}.dt_0_4.txt"), sep="\t")

    # --- estimability: Stage B contrasts perturb vs NTC; needs non-zero usage variance among NTC ---
    meta = pd.read_csv(args.metadata, index_col=0)
    ntc = meta.index[meta["gene"] == "non-targeting"]
    ntc = [c for c in ntc if c in U_df.index]
    est = pd.DataFrame({
        "program": [int(c) for c in U_df.columns],
        "nz_all": [(U_df[c] > 0).sum() for c in U_df.columns],
        "nz_ctrl": [(U_df.loc[ntc, c] > 0).sum() for c in U_df.columns],
        "sd_ctrl": [float(U_df.loc[ntc, c].std()) for c in U_df.columns],
        "rel_err_dummy": np.nan,
    }).set_index("program")
    est["estimable"] = est.sd_ctrl > 0
    est.to_csv(os.path.join(args.results, "cNMF", args.ct, "program_estimability.tsv"), sep="\t")
    alive = [str(p) for p in est.index[est.estimable]]
    with open(os.path.join(args.results, "cNMF", args.ct, "alive_programs.txt"), "w") as f:
        f.write("\n".join(alive) + "\n")
    print(f"[proj] {args.ct}: estimable {len(alive)}/{K}; NOT estimable {list(est.index[~est.estimable])}"
          f" -> {out}", flush=True)


if __name__ == "__main__":
    main()
