#!/usr/bin/env python3
"""
06_prepare_inputs_backman.py
----------------------------
Build the GeneBayes inputs for the Backman-2021 RA panel:
  * traits/<label>.summary_statistics.csv   (ensg,beta,standard_error)  -- response
  * genelists/<label>.all.txt               (ENSG in response AND features)
Reuses the shared numeric feature table built by 02_prepare_inputs_panel.py.

Backman beta = ln(odds_ratio) (log-odds); standard_error is reported directly and was
verified consistent with the reported p (no reconstruction). SE cap 100 is inert here
(max observed SE ~19.5) but kept for parity with the Genebass path.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "LoF" / "raw_burden"
FEAT_OUT = ROOT / "features" / "gene_features_ensg_numeric.tsv"
TRAIT_DIR = ROOT / "traits"; TRAIT_DIR.mkdir(exist_ok=True)
GL_DIR = ROOT / "genelists"; GL_DIR.mkdir(exist_ok=True)

PANEL = [
    "Backman_2021_RA",
    "Backman_2021_RA_M06",
    "Backman_2021_RA_alt",
    "Backman_2021_RA_M05_seropos",
    "Backman_2021_RA_M060_seroneg",
]
SE_CAP = 100.0

if not FEAT_OUT.exists():
    raise SystemExit(f"missing {FEAT_OUT}; run 02_prepare_inputs_panel.py first")
with open(FEAT_OUT) as f:
    next(f)
    feat_genes = {ln.split("\t", 1)[0] for ln in f}
print(f"[features] {len(feat_genes)} genes in shared feature table")

for label in PANEL:
    raw = RAW_DIR / f"{label}_M1_001.summary_statistics.csv"
    resp = {}
    for r in csv.DictReader(open(raw)):
        resp[r["ensg"]] = (float(r["beta"]), min(float(r["standard_error"]), SE_CAP))
    out = TRAIT_DIR / f"{label}.summary_statistics.csv"
    with open(out, "w") as fh:
        fh.write("ensg,beta,standard_error\n")
        for ensg, (b, se) in resp.items():
            fh.write(f"{ensg},{b:.10g},{se:.10g}\n")
    common = sorted(set(resp) & feat_genes)
    (GL_DIR / f"{label}.all.txt").write_text("\n".join(common) + "\n")
    print(f"[{label}] response {len(resp)} -> {out.name}; response∩features {len(common)}")
