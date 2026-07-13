#!/usr/bin/env python3
"""
02_prepare_inputs_panel.py
--------------------------
For each phenotype in the RA panel, build the GeneBayes inputs the paper's
genebayes.lof.gene_level.py reads:
  * traits/<name>.summary_statistics.csv   (ensg,beta,standard_error)  -- response
  * genelists/<name>.all.txt               (ENSG present in response AND features)
The shared numeric feature table (features/gene_features_ensg_numeric.tsv,
ensg + 1247 numeric cols, hgnc/chrom dropped) is built once from the GeneBayes
curated feature file and reused across all phenotypes.
"""
import gzip
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "LoF" / "raw_burden"
FEAT_GZ = ROOT / "features" / "gene_features_for_s_het.tsv.gz"
FEAT_OUT = ROOT / "features" / "gene_features_ensg_numeric.tsv"
TRAIT_DIR = ROOT / "traits"; TRAIT_DIR.mkdir(exist_ok=True)
GL_DIR = ROOT / "genelists"; GL_DIR.mkdir(exist_ok=True)

PANEL = ["Genebass_RF", "Genebass_RA_custom", "Genebass_RA_M06"]
SE_CAP = 100.0


def build_feature_table():
    """ensg + numeric features (drop hgnc, chrom). Returns set of feature ENSGs."""
    feat_genes = set()
    with gzip.open(FEAT_GZ, "rt") as fin, open(FEAT_OUT, "w") as fout:
        hdr = fin.readline().rstrip("\n").split("\t")
        drop = {"hgnc", "chrom"}
        keep_idx = [i for i, c in enumerate(hdr) if c not in drop]
        ensg_col = hdr.index("ensg")
        new_hdr = ["ensg"] + [hdr[i] for i in keep_idx if hdr[i] != "ensg"]
        fout.write("\t".join(new_hdr) + "\n")
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            feat_genes.add(parts[ensg_col])
            vals = [parts[i] for i in keep_idx if hdr[i] != "ensg"]
            fout.write(parts[ensg_col] + "\t" + "\t".join(vals) + "\n")
    print(f"[features] {len(feat_genes)} genes x {len(new_hdr)-1} numeric features -> {FEAT_OUT.name}")
    return feat_genes


def prepare_trait(name, feat_genes):
    raw = RAW_DIR / f"{name}_M1_pLoF.summary_statistics.csv"
    resp = {}
    with open(raw) as fh:
        for r in csv.DictReader(fh):
            resp[r["ensg"]] = (float(r["beta"]), min(float(r["standard_error"]), SE_CAP))
    resp_out = TRAIT_DIR / f"{name}.summary_statistics.csv"
    with open(resp_out, "w") as fh:
        fh.write("ensg,beta,standard_error\n")
        for ensg, (beta, se) in resp.items():
            fh.write(f"{ensg},{beta:.10g},{se:.10g}\n")
    common = sorted(set(resp) & feat_genes)
    gl_out = GL_DIR / f"{name}.all.txt"
    gl_out.write_text("\n".join(common) + "\n")
    print(f"[{name}] response {len(resp)} genes -> {resp_out.name}; "
          f"response∩features {len(common)} -> {gl_out.name}")


if __name__ == "__main__":
    if not FEAT_OUT.exists():
        feat_genes = build_feature_table()
    else:
        with open(FEAT_OUT) as f:
            next(f)
            feat_genes = {ln.split("\t", 1)[0] for ln in f}
        print(f"[features] reuse existing {FEAT_OUT.name} ({len(feat_genes)} genes)")
    for name in PANEL:
        prepare_trait(name, feat_genes)
