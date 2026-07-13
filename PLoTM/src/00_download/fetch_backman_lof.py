#!/usr/bin/env python3
"""
05_fetch_backman_RA.py
----------------------
Fetch the rheumatoid-arthritis **gene-based burden** summary statistics from
**Backman et al. 2021** (UK Biobank 454,787 exomes) — the paper's own LoF source —
via the GWAS Catalog, and write them in the exact schema/naming the Ota et al. (2026)
pipeline uses:  data/LoF/raw_burden/Backman_2021_<label>_M1_001.summary_statistics.csv

Provenance: Backman 2021's Data Availability states "Summary statistics for the rare
variants tested in this study are also available in the GWAS Catalog ... listed
separately for single variants and burden tests." Burden studies are the ones whose
trait name ends in "(Gene-based burden)".

Mask: per the GWAS Catalog meta.yaml for these studies —
    M1 = putative loss-of-function (pLoF) only;  M3 = pLoF + deleterious missense
    AAF rollups:  M1.singleton | M1.0001 (MAF<=0.001%) | M1.001 (MAF<=0.01%)
                  M1.01 (MAF<=0.1%) | M1.1 (MAF<=1%)
The paper's files are named `Backman_2021_<id>_M1_001` => mask **M1.001** (pLoF, MAF<=0.01%).
We select exactly that mask.

Effect: these RA traits are binary, so the file reports `odds_ratio`; the burden
effect on the log-odds scale is beta = ln(OR), and `standard_error` is already the
SE of that log-OR (regenie ADD-WGR-FIRTH). **No SE reconstruction is needed** —
unlike the Genebass route.
"""
import csv
import gzip
import io
import math
import subprocess
import sys
from pathlib import Path

MASK = "M1.001"  # pLoF, MAF <= 0.01%  (the paper's _M1_001)

# label,                              accession,      description
PANEL = [
    ("Backman_2021_RA",               "GCST90085490", "Rheumatoid arthritis (curated), 8,561 cases"),
    ("Backman_2021_RA_M06",           "GCST90084378", "ICD10 M06 Other rheumatoid arthritis, 6,100 cases"),
    ("Backman_2021_RA_alt",           "GCST90081859", "Rheumatoid arthritis (alt. definition), 5,095 cases"),
    ("Backman_2021_RA_M05_seropos",   "GCST90084373", "ICD10 M05 RA with rheumatoid factor (seropositive), 769 cases"),
    ("Backman_2021_RA_M060_seroneg",  "GCST90084374", "ICD10 M06.0 RA without rheumatoid factor (seronegative), 706 cases"),
]

ROOT = Path(__file__).resolve().parent.parent
RAW_DIR = ROOT / "data" / "LoF" / "raw_burden"
RAW_DIR.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / "data" / "LoF" / "backman_cache"
CACHE.mkdir(parents=True, exist_ok=True)


def ftp_url(acc: str) -> str:
    n = int(acc.replace("GCST", ""))
    lo = (n - 1) // 1000 * 1000 + 1
    hi = lo + 999
    return (f"https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            f"GCST{lo:08d}-GCST{hi:08d}/{acc}/{acc}_buildGRCh38.tsv.gz")


def norm_sf(z: float) -> float:
    """two-sided p from |z|"""
    return math.erfc(abs(z) / math.sqrt(2.0))


def fetch(acc: str) -> Path:
    dest = CACHE / f"{acc}.tsv.gz"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    subprocess.run(["curl", "-sL", "--retry", "5", "--retry-delay", "5",
                    "--max-time", "600", "-o", str(dest), ftp_url(acc)], check=True)
    return dest


def process(label, acc, desc):
    src = fetch(acc)
    rows = []
    n_mask = 0
    consist = []
    with gzip.open(src, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            name = r["Name"]
            if not name.endswith(f".GENE.{MASK}"):
                continue
            n_mask += 1
            # Name = SYMBOL(ENSG).GENE.M1.001
            head = name.split(".GENE.")[0]
            if "(" not in head or not head.endswith(")"):
                continue
            sym, ensg = head[:-1].split("(", 1)
            try:
                orr = float(r["odds_ratio"]); se = float(r["standard_error"])
            except (ValueError, TypeError, KeyError):
                continue
            if not (orr > 0) or not (se > 0) or not math.isfinite(orr) or not math.isfinite(se):
                continue
            beta = math.log(orr)          # burden effect on log-odds scale
            z = beta / se
            rows.append((ensg, sym, f"chr{r['chromosome']}", beta, se, z))
            # consistency: does reported p match Wald p from beta/se?
            try:
                p = float(r["p_value"])
                if 1e-12 < p < 1:
                    consist.append((norm_sf(z), p))
            except (ValueError, TypeError, KeyError):
                pass

    rows.sort(key=lambda x: abs(x[5]), reverse=True)
    out = RAW_DIR / f"{label}_M1_001.summary_statistics.csv"
    with open(out, "w") as fh:
        fh.write("ensg,gene,chrom,beta,standard_error,z\n")
        for ensg, sym, chrom, beta, se, z in rows:
            fh.write(f"{ensg},{sym},{chrom},{beta:.10g},{se:.10g},{z:.10g}\n")

    maxz = max((abs(r[5]) for r in rows), default=0.0)
    # median |log10 p_wald / p_reported| as a sanity metric
    if consist:
        ratios = [abs(math.log10(a / b)) for a, b in consist if a > 0 and b > 0]
        ratios.sort()
        med = ratios[len(ratios) // 2]
    else:
        med = float("nan")
    print(f"[{label}] {desc}")
    print(f"    acc={acc}  mask={MASK}  rows_with_mask={n_mask}  usable_genes={len(rows)}")
    print(f"    max|Z|={maxz:.2f}  #|Z|>3={sum(abs(r[5])>3 for r in rows)}  "
          f"median|log10(p_wald/p_reported)|={med:.3f}")
    print(f"    -> {out.name}")


if __name__ == "__main__":
    for label, acc, desc in PANEL:
        process(label, acc, desc)
