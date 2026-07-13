#!/usr/bin/env python3
"""
12_qc_converged.py — final QC over the CONVERGED (`*_conv`) posteriors for all 8 RA traits,
with a 50-iter-vs-converged comparison and the seronegative degeneracy test.

Reports per trait: gamma distribution, shrinkage, prior<0, sign-flip vs beta (overall and
among informative |Z|>2 genes), corr(gamma_50, gamma_conv), and early-stopping best iter.
Then: within-Backman gamma agreement, cross-source concordance, canonical RA genes.
"""
import csv, math, re, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
POST = ROOT / "data" / "LoF" / "GeneBayes_posterior"
RAW = ROOT / "data" / "LoF" / "raw_burden"
LOGS = ROOT / "genebayes_run"

# label: (posterior stem, raw burden file)
TRAITS = {
    "GB_RF":          ("Genebass_RF",                  "Genebass_RF_M1_pLoF"),
    "GB_RA_custom":   ("Genebass_RA_custom",           "Genebass_RA_custom_M1_pLoF"),
    "GB_RA_M06":      ("Genebass_RA_M06",              "Genebass_RA_M06_M1_pLoF"),
    "BK_RA_curated":  ("Backman_2021_RA",              "Backman_2021_RA_M1_001"),
    "BK_RA_M06":      ("Backman_2021_RA_M06",          "Backman_2021_RA_M06_M1_001"),
    "BK_RA_alt":      ("Backman_2021_RA_alt",          "Backman_2021_RA_alt_M1_001"),
    "BK_seropos":     ("Backman_2021_RA_M05_seropos",  "Backman_2021_RA_M05_seropos_M1_001"),
    "BK_seroneg":     ("Backman_2021_RA_M060_seroneg", "Backman_2021_RA_M060_seroneg_M1_001"),
}
CANON = ["TNFAIP3", "IL6R", "PTPN22", "TYK2", "STAT4", "CTLA4", "IRF5", "REL", "NFKBIA"]


def pear(x, y):
    n = len(x)
    if n < 3: return float("nan")
    mx, my = sum(x)/n, sum(y)/n
    c = sum((a-mx)*(b-my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a-mx)**2 for a in x)); sy = math.sqrt(sum((b-my)**2 for b in y))
    return c/(sx*sy) if sx and sy else float("nan")


def load_post(stem):
    p = POST / f"{stem}.per_gene_estimates.tsv"
    if not p.exists(): return None, None
    g, pr = {}, {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["post_mean"] in ("", "nan"): continue
        g[r["ensg"]] = float(r["post_mean"]); pr[r["ensg"]] = float(r["prior_mean"])
    return g, pr


def best_iter(stem):
    lg = LOGS / f"{stem}.run.log"
    if not lg.exists(): return "-"
    m = None
    for line in open(lg, errors="replace"):
        mm = re.search(r"Best iteration / VAL(\d+)", line)
        if mm: m = mm.group(1)
    return m or "no-early-stop"


gam, pri, beta, zz, sym = {}, {}, {}, {}, {}
flags = {}
print(f"{'trait':<15}{'n':>7}{'gam<0':>8}{'prior<0':>9}{'shrunk':>8}{'flip':>7}{'flip|Z|>2':>11}"
      f"{'r(g50,gconv)':>14}{'bestIter':>10}")
for lab, (stem, rawf) in TRAITS.items():
    g, pr = load_post(stem + "_conv")
    if g is None:
        print(f"{lab:<15}  (converged posterior missing — skipping)"); continue
    g50, _ = load_post(stem)                      # 50-iter version, may be absent
    b, z, s = {}, {}, {}
    for r in csv.DictReader(open(RAW / f"{rawf}.summary_statistics.csv")):
        b[r["ensg"]] = float(r["beta"]); z[r["ensg"]] = float(r["z"]); s[r["ensg"]] = r["gene"]
    gam[lab], pri[lab], beta[lab], zz[lab] = g, pr, b, z; sym.update(s)

    v = list(g.values()); pv = list(pr.values())
    sh = [e for e in g if e in b and b[e] != 0]
    flip = sum((g[e] > 0) != (b[e] > 0) for e in sh) / len(sh)
    strong = [e for e in sh if abs(z[e]) > 2]
    flips = (sum((g[e] > 0) != (b[e] > 0) for e in strong) / len(strong)) if strong else float("nan")
    shr = sum(abs(g[e]) < abs(b[e]) for e in sh) / len(sh)
    if g50:
        c = sorted(set(g) & set(g50)); r50 = pear([g[e] for e in c], [g50[e] for e in c])
    else:
        r50 = float("nan")
    print(f"{lab:<15}{len(v):>7}{sum(x<0 for x in v)/len(v):>7.1%}{sum(x<0 for x in pv)/len(pv):>8.1%}"
          f"{shr:>8.0%}{flip:>7.1%}{flips:>11.1%}{r50:>14.3f}{best_iter(stem+'_conv'):>10}")

# ---- prior-domination / sign-degeneracy screen + cross-source scale --------------
# When SE >> |beta| the likelihood is flat, the posterior collapses onto the feature-based
# prior, and sign(gamma) carries no trait-specific direction. Screen BEFORE any sign test.
print("\n=== prior-domination screen (run before any directional/sign test) ===")
print(f"{'trait':<15}{'sd(gam)':>10}{'p99|gam|':>10}{'gam<0':>8}{'r(post,prior)':>15}{'medSE/med|b|':>14}  flag")
for lab in gam:
    g = gam[lab]; pr = pri[lab]
    ens = list(g)
    gv = [g[e] for e in ens]; pv = [pr[e] for e in ens]
    rpp = pear(gv, pv)
    fneg = sum(x < 0 for x in gv)/len(gv)
    rawf = TRAITS[lab][1]
    ab, se = [], []
    for r in csv.DictReader(open(RAW / f"{rawf}.summary_statistics.csv")):
        ab.append(abs(float(r["beta"]))); se.append(float(r["standard_error"]))
    ratio = st.median(se)/st.median(ab)
    p99 = sorted(abs(x) for x in gv)[int(0.99*len(gv))]
    fl = []
    if rpp > 0.95: fl.append("PRIOR-DOMINATED")
    if fneg < 0.02 or fneg > 0.98: fl.append("SIGN-DEGENERATE")
    flags[lab] = fl
    print(f"{lab:<15}{st.pstdev(gv):>10.4f}{p99:>10.4f}{fneg:>7.1%}{rpp:>15.3f}{ratio:>14.2f}  {', '.join(fl) or 'ok'}")
bad = [l for l, f in flags.items() if f]
if bad:
    print(f"  -> DO NOT use for sign-based tests: {', '.join(bad)}")
print("  NOTE: gamma scale differs ~30-60x between Backman (ln OR) and Genebass (burden score).")
print("        Never share a fixed |gamma| cutoff across sources; use a percentile-matched one.")

bk = [l for l in gam if l.startswith("BK_")]
if len(bk) >= 2:
    print("\n=== within-Backman converged gamma correlation ===")
    print(" "*15 + "".join(f"{l.replace('BK_',''):>12}" for l in bk))
    for a in bk:
        row = "".join(f"{pear([gam[a][e] for e in sorted(set(gam[a])&set(gam[b]))],[gam[b][e] for e in sorted(set(gam[a])&set(gam[b]))]):>12.3f}" for b in bk)
        print(f"  {a.replace('BK_',''):<13}{row}")

print("\n=== raw-Z vs converged-gamma agreement (same pairs) ===")
for a, b in [("BK_RA_curated","BK_RA_M06"), ("BK_seropos","BK_seroneg"),
             ("BK_RA_curated","GB_RA_custom"), ("BK_RA_M06","GB_RA_M06")]:
    if a in gam and b in gam:
        sz = sorted(set(zz[a]) & set(zz[b])); sg = sorted(set(gam[a]) & set(gam[b]))
        print(f"  {a:<15} vs {b:<15} r(Z)={pear([zz[a][e] for e in sz],[zz[b][e] for e in sz]):+.3f}"
              f"   r(gamma)={pear([gam[a][e] for e in sg],[gam[b][e] for e in sg]):+.3f}")

if gam:
    labs = list(gam)
    print("\n=== canonical RA genes: converged gamma ===")
    s2e = {v: k for k, v in sym.items()}
    print(f"  {'gene':<9}" + "".join(f"{l:>14}" for l in labs))
    for gname in CANON:
        e = s2e.get(gname)
        if not e: continue
        print(f"  {gname:<9}" + "".join(f"{gam[l].get(e, float('nan')):>14.3f}" for l in labs))
print("\nDone.")
