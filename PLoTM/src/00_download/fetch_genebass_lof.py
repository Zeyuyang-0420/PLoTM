#!/usr/bin/env python3
"""
01_fetch_genebass_panel.py
--------------------------
Fetch predicted-loss-of-function (pLoF) *gene burden* association statistics for a
curated panel of rheumatoid-arthritis phenotypes from Genebass (UK Biobank,
394,841 exomes; Karczewski et al. 2022) and write each in the Ota et al. (2026)
pipeline's UKB-LoF schema (same layout as the paper's Backman `*_M1_001.summary_statistics.csv`):

    ensg, gene, chrom, beta, standard_error, z

Panel (complementary facets of RA):
  * Genebass_RF        Rheumatoid factor (irnt)         33,546  continuous  <- well-powered endophenotype anchor
  * Genebass_RA_custom Rheumatoid arthritis (curated)    6,101  case/control <- primary RA diagnosis
  * Genebass_RA_M06    M06 other RA (ICD first-occur.)   7,852  case/control <- independent RA diagnosis (replication)

SE reconstruction (public API exposes BETA_Burden + Pvalue_Burden, not SE):
    |Z| = Phi^{-1}(1 - P/2)  ;  SE = |beta| / |Z|
the exact Wald relation the paper's own scripts use in reverse
(Figure2/2_volcanoplot.R, Figure3/2_GWAS_enrich.R). Cleanest for the continuous
RF trait (linear model, no saddlepoint approximation).
"""
import json
import math
import subprocess
import sys
from pathlib import Path
from statistics import NormalDist

API = "https://main.genebass.org/api"
BURDEN_SET = "pLoF"   # predicted LoF (== Backman M1 mask)

PANEL = [
    # out_name,           analysis_id,                                   description
    ("Genebass_RF",        "continuous-30820-both_sexes--irnt",          "Rheumatoid factor (irnt), 33,546"),
    ("Genebass_RA_custom", "categorical-RA_custom-both_sexes--custom",   "Rheumatoid arthritis (RA), curated, 6,101 cases"),
    ("Genebass_RA_M06",    "icd_first_occurrence-131850-both_sexes--",   "M06 other rheumatoid arthritis (ICD), 7,852 cases"),
]

HERE = Path(__file__).resolve().parent
LOF = HERE.parent / "data" / "LoF"
RAW_DIR = LOF / "raw_burden"; RAW_DIR.mkdir(parents=True, exist_ok=True)
VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y"}


def curl_json(url):
    res = subprocess.run(
        ["curl", "-s", "--compressed", "--retry", "5", "--retry-delay", "5",
         "--max-time", "180", url],
        capture_output=True)
    if res.returncode != 0:
        sys.exit(f"curl failed for {url}: {res.stderr.decode()[:300]}")
    return json.loads(res.stdout)


def fetch_one(out_name, analysis_id, desc):
    url = f"{API}/analysis/{analysis_id}/gene-manhattan?burdenSet={BURDEN_SET}"
    genes = curl_json(url)
    (LOF / f"{out_name}_pLoF.gene-manhattan.json").write_text(json.dumps(genes))
    rows = []
    n_bad = 0
    for g in genes:
        ensg, beta, pb = g.get("gene_id"), g.get("BETA_Burden"), g.get("Pvalue_Burden")
        chrom, sym = g.get("chrom"), g.get("gene_symbol")
        if ensg is None or beta is None or pb is None:
            n_bad += 1
            continue
        pb = min(max(float(pb), 1e-300), 1.0 - 1e-15)
        z_abs = -NormalDist().inv_cdf(pb / 2.0)          # Phi^{-1}(1 - P/2)
        if z_abs <= 0 or not math.isfinite(z_abs):
            n_bad += 1
            continue
        se = abs(beta) / z_abs
        z_signed = beta / se if se > 0 else 0.0
        chrom_fmt = f"chr{chrom}" if str(chrom) in VALID_CHROMS else str(chrom)
        rows.append((ensg, sym, chrom_fmt, beta, se, z_signed))
    rows.sort(key=lambda r: abs(r[5]), reverse=True)
    out = RAW_DIR / f"{out_name}_M1_pLoF.summary_statistics.csv"
    with open(out, "w") as fh:
        fh.write("ensg,gene,chrom,beta,standard_error,z\n")
        for ensg, sym, chrom, beta, se, z in rows:
            fh.write(f"{ensg},{sym},{chrom},{beta:.10g},{se:.10g},{z:.10g}\n")
    maxz = max((abs(r[5]) for r in rows), default=0.0)
    print(f"[{out_name}] {desc}")
    print(f"    {len(rows)} genes (dropped {n_bad}); max|Z|={maxz:.2f}; "
          f"#|Z|>3={sum(abs(r[5])>3 for r in rows)} -> {out.name}")
    return rows


if __name__ == "__main__":
    for name, aid, desc in PANEL:
        fetch_one(name, aid, desc)
