#!/usr/bin/env python
"""build_demo_data.py — flatten pipeline_results/ into PLoTM_explorer/data/demo_data.json and
copy the display figures into static/figures/<disease>/.

The output schema is DISEASE-KEYED so a second disease (e.g. SLE) drops in with a new entry in
DISEASES below and a re-run — no code change, no UI change. RA is the reference/flagship.

Run with the cNMF env's python (needs pandas); always with PYTHONNOUSERSITE=1:
    PYTHONNOUSERSITE=1 /mnt/scratch/ZY2/.envs/perturbseq-repro/bin/python build_demo_data.py --disease RA
"""
import argparse, glob, json, math, os, re, shutil, sys
from datetime import datetime, timezone

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("PROGLOF_ROOT", os.path.abspath(os.path.join(HERE, "..")))  # .../Hackathon

# ---------------------------------------------------------------------------
# Per-disease source map. Adding a disease = adding an entry here (SLE-ready).
# ---------------------------------------------------------------------------
DISEASES = {
    "RA": {
        "name": "Rheumatoid arthritis",
        "cell_type": "CD4⁺ T cells — genome-wide CRISPRi Perturb-seq (GWCD4i)",
        "results_root": os.path.join(ROOT, "pipeline_results"),
        "spectra": {
            "rest": os.path.join(ROOT, "tcell_perturbseq/cNMF_RA_analysis/cNMF/"
                    "GWCD4i_rest_pseudobulk/test1/test1.gene_spectra_score.k_60.dt_0_4.txt"),
            "stim48hr": os.path.join(ROOT, "tcell_perturbseq/cNMF_RA_analysis/cNMF/"
                    "GWCD4i_stim48_pseudobulk/test1/test1.gene_spectra_score.k_60.dt_0_4.txt"),
        },
        "trait_desc": {
            "Genebass_RF_conv": "rheumatoid factor (continuous endophenotype)",
            "Genebass_RA_custom_conv": "RA, custom case definition",
            "Genebass_RA_M06_conv": "RA, phecode M06",
            "Backman_2021_RA_conv": "RA",
            "Backman_2021_RA_alt_conv": "RA, alternate definition",
            "Backman_2021_RA_M06_conv": "RA, phecode M06",
            "Backman_2021_RA_M05_seropos_conv": "seropositive RA (M05)",
            "Backman_2021_RA_M060_seroneg_conv": "seronegative RA (M060)",
        },
        # trait groups for the combined interactive network (mirror combined_figures_annotated;
        # the 2 QC-degenerate Backman facets are intentionally excluded).
        "groups": {
            "Group1_RA": {"label": "Group 1 · RA (5 phenotypes)",
                          "traits": ["Backman_2021_RA_conv", "Backman_2021_RA_M06_conv",
                                     "Backman_2021_RA_alt_conv", "Genebass_RA_custom_conv",
                                     "Genebass_RA_M06_conv"]},
            "Group2_RF": {"label": "Group 2 · RF (endophenotype)",
                          "traits": ["Genebass_RF_conv"]},
        },
    },
}

# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def canon(cond):
    """tidy uses 'stim48'; summary/figures use 'stim48hr'. Normalise to {rest, stim48hr}."""
    return "stim48hr" if str(cond).lower().startswith("stim") else "rest"

def fnum(x):
    """float or None (NaN/empty/unparseable -> None) so the JSON is valid (no NaN literals)."""
    try:
        v = float(x)
        return None if math.isnan(v) or math.isinf(v) else v
    except (TypeError, ValueError):
        return None

def sstr(x):
    if x is None:
        return ""
    s = str(x)
    return "" if s.lower() in ("nan", "none") else s

def read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)

def find_gencode(results_root):
    hits = glob.glob(os.path.join(results_root, "**", "gencode_v41_gname_gid_ALL_sorted_onlyID"),
                     recursive=True)
    return hits[0] if hits else None

def load_gencode(path):
    m = {}
    if not path or not os.path.exists(path):
        return m
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                m[p[0]] = p[1]
    return m

def top_genes_per_program(spectra_path, ensg2sym, topn=12):
    """gene_spectra_score: rows = programs 1..60, cols = ENSG. Return {prog_int: [symbols]}."""
    out = {}
    if not spectra_path or not os.path.exists(spectra_path):
        return out
    df = pd.read_csv(spectra_path, sep="\t", index_col=0)
    for prog, row in df.iterrows():
        try:
            pi = int(prog)
        except (TypeError, ValueError):
            continue
        top = row.sort_values(ascending=False).head(topn)
        out[pi] = [ensg2sym.get(g, g) for g in top.index]
    return out

def prog_num(p):
    m = re.match(r"P(\d+)", str(p))
    return int(m.group(1)) if m else 9999


def _iint(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


def _glob1(pattern):
    hits = glob.glob(pattern)
    return hits[0] if hits else None


def _read_map(rr, trait, cond):
    """Fig-5a graph for one (trait, condition): (nodes, edges) DataFrames. Missing/empty -> empty."""
    def _rd(path, cols):
        if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
            return pd.DataFrame(columns=cols)
        try:
            return pd.read_csv(path, sep="\t", dtype=str)
        except pd.errors.EmptyDataError:
            return pd.DataFrame(columns=cols)
    nd = _glob1(os.path.join(rr, trait, "tables", f"{cond}__*_map_nodes.tsv"))
    ed = _glob1(os.path.join(rr, trait, "tables", f"{cond}__*_map_edges.tsv"))
    return _rd(nd, ["id", "kind", "sign"]), _rd(ed, ["GENE", "PROGRAM", "beta", "padj", "type"])


def _selected(rr, trait, cond):
    f = _glob1(os.path.join(rr, trait, "tables", f"{cond}__SelectedPrograms_*.txt"))
    if not f:
        return [], []
    df = pd.read_csv(f, sep="\t", dtype=str)
    if not len(df):
        return [], []
    ps = [p for p in str(df.iloc[0].get("Program_selected", "")).split(",") if p]
    rs = [p for p in str(df.iloc[0].get("Regulator_selected", "")).split(",") if p]
    return ps, rs


def build_networks(rr, cfg, ann_by, meta_sig, progMeta_by, short_by):
    """Aggregate the per-trait Fig-5a tri-partite graphs into one interactive graph per
    (condition, trait-group): genes -> programs -> traits. Mirrors combined_figures_annotated."""
    groups = cfg.get("groups", {})
    out = {"rest": {}, "stim48hr": {}}
    for cond in ("rest", "stim48hr"):
        for gid, g in groups.items():
            gene_sign, gene_traits, prog_traits, prog_role, edge_agg = {}, {}, {}, {}, {}
            trait_ids = []
            for trait in g["traits"]:
                short = short_by.get(trait, trait)
                trait_ids.append(short)
                nodes, edges = _read_map(rr, trait, cond)
                for _, r in nodes.iterrows():
                    if r["kind"] == "gene":
                        gene_sign.setdefault(r["id"], []).append(_iint(r.get("sign")))
                        gene_traits.setdefault(r["id"], set()).add(short)
                    elif r["kind"] == "program":
                        prog_traits.setdefault(r["id"], set()).add(short)
                ps, rs = _selected(rr, trait, cond)
                for p in ps:
                    prog_role.setdefault(p, set()).add("program"); prog_traits.setdefault(p, set()).add(short)
                for p in rs:
                    prog_role.setdefault(p, set()).add("regulator"); prog_traits.setdefault(p, set()).add(short)
                for _, r in edges.iterrows():
                    key = (r["GENE"], r["PROGRAM"])
                    e = edge_agg.setdefault(key, {"type": r.get("type"), "betas": [], "traits": set()})
                    b = fnum(r.get("beta"))
                    if b is not None:
                        e["betas"].append(b)
                    e["traits"].add(short)
                    if r.get("type") == "regulates":
                        e["type"] = "regulates"
                    gene_traits.setdefault(r["GENE"], set()).add(short)
            gene_is_reg = {gg for (gg, _), e in edge_agg.items() if e["type"] == "regulates"}
            nodes_out, edges_out = [], []
            for gene, traits in gene_traits.items():
                signs = [s for s in gene_sign.get(gene, []) if s is not None]
                sign = max(set(signs), key=signs.count) if signs else 0
                nodes_out.append({"id": gene, "type": "gene",
                                  "role": "regulator" if gene in gene_is_reg else "member",
                                  "sign": sign, "n_traits": len(traits), "traits": sorted(traits)})
            for prog, traits in prog_traits.items():
                roles = prog_role.get(prog, set())
                role = "both" if len(roles) == 2 else (next(iter(roles)) if roles else "program")
                ann = ann_by.get(cond, {}).get(prog, "")
                nodes_out.append({"id": prog, "type": "program", "annotation": ann,
                                  "annotated": bool(ann), "role": role, "n_traits": len(traits),
                                  "traits": sorted(traits), "sig": prog in meta_sig.get(cond, set())})
            for tn in sorted(set(trait_ids)):
                nodes_out.append({"id": tn, "type": "trait"})
            for (gene, prog), e in edge_agg.items():
                beta = max(e["betas"], key=abs) if e["betas"] else None
                edges_out.append({"source": gene, "target": prog, "cls": e["type"],
                                  "beta": beta, "n_traits": len(e["traits"]), "traits": sorted(e["traits"])})
            for prog, traits in prog_traits.items():
                for short in traits:
                    meta = progMeta_by.get((cond, prog, short)) or {}
                    w = meta.get("P")
                    d = meta.get("dir")
                    edges_out.append({"source": prog, "target": short, "cls": "selects",
                                      "weight": (-math.log10(w) if (w and w > 0) else 1.0),
                                      "dir": (1 if (d or 0) > 0 else -1),
                                      "effect": meta.get("effect")})
            out[cond][gid] = {"label": g["label"], "nodes": nodes_out, "edges": edges_out,
                              "n_programs": sum(1 for n in nodes_out if n["type"] == "program"),
                              "n_genes": sum(1 for n in nodes_out if n["type"] == "gene")}
    return out


def build_regulators(rr, all_traits, ann_by, short_by, flag_by_short):
    """One entry per regulator gene (a gene with a `regulates` edge anywhere): which programs it
    regulates (β, direction, annotation), across which traits/conditions. Surfaces shared regulators."""
    reg, gene_sign = {}, {}
    for trait in all_traits:
        short = short_by.get(trait, trait)
        flag = flag_by_short.get(short, "ok")
        for cond in ("rest", "stim48hr"):
            nodes, edges = _read_map(rr, trait, cond)
            for _, r in nodes.iterrows():
                if r["kind"] == "gene":
                    gene_sign.setdefault(r["id"], []).append(_iint(r.get("sign")))
            for _, r in edges.iterrows():
                if r.get("type") != "regulates":
                    continue
                gene, prog, b = r["GENE"], r["PROGRAM"], fnum(r.get("beta"))
                g = reg.setdefault(gene, {"targets": [], "traits": set(), "conditions": set(), "programs": set()})
                g["targets"].append({"program": prog, "annotation": ann_by.get(cond, {}).get(prog, ""),
                                     "beta": b, "dir": (1 if (b or 0) > 0 else -1),
                                     "trait": short, "condition": cond, "flag": flag})
                g["traits"].add(short); g["conditions"].add(cond); g["programs"].add(prog)
    out = []
    for gene, g in reg.items():
        signs = [s for s in gene_sign.get(gene, []) if s is not None]
        sign = max(set(signs), key=signs.count) if signs else 0
        clean_traits = sorted({t["trait"] for t in g["targets"] if t["flag"] == "ok"})
        out.append({"gene": gene, "sign": sign, "n_programs": len(g["programs"]),
                    "n_targets": len(g["targets"]), "n_traits": len(g["traits"]),
                    "traits": sorted(g["traits"]), "conditions": sorted(g["conditions"]),
                    "n_clean_traits": len(clean_traits), "clean_traits": clean_traits,
                    "shared": len(clean_traits) >= 2,   # shared = recurs across QC-clean phenotypes
                    "targets": sorted(g["targets"], key=lambda t: -abs(t["beta"] or 0))})
    # rank by clean-phenotype recurrence first, so real cross-trait regulators lead (not noise-program
    # hits that only recur among QC-flagged phenotypes).
    out.sort(key=lambda r: (-r["n_clean_traits"], -r["n_traits"], -r["n_programs"], r["gene"]))
    return out

# ---------------------------------------------------------------------------
# build one disease
# ---------------------------------------------------------------------------
def build_disease(did, cfg):
    rr = cfg["results_root"]
    tidy = read_tsv(os.path.join(rr, "condition_comparison/tables/tidy_results.tsv"))
    meta = read_tsv(os.path.join(rr, "condition_comparison/tables/meta_regulator_RA5.tsv"))
    thr = read_tsv(os.path.join(rr, "lof_thresholds.tsv"))
    summ = read_tsv(os.path.join(rr, "summary.tsv"))

    tidy["cond"] = tidy["condition"].map(canon)
    meta["cond"] = meta["condition"].map(canon)
    summ["cond"] = summ["cond"].map(canon)

    ensg2sym = load_gencode(find_gencode(rr))
    topg = {c: top_genes_per_program(cfg["spectra"].get(c), ensg2sym) for c in ("rest", "stim48hr")}

    thr_by = {r["trait"]: r for _, r in thr.iterrows()}

    # ---- traits[] (metadata + per-condition panel headline from summary.tsv) ----
    traits = []
    trait_order = list(cfg["trait_desc"].keys())
    seen = set(tidy["trait"])
    for tid in trait_order:
        if tid not in seen:
            continue
        sub = tidy[tidy["trait"] == tid].iloc[0]
        tr = thr_by.get(tid, {})
        panels = {}
        for _, s in summ[summ["trait"] == tid].iterrows():
            panels[s["cond"]] = {
                "fig4c_prog_hit": sstr(s.get("prog_bonf")),
                "fig4c_reg_hit": sstr(s.get("reg_bonf")),
                "fig5a_perm_P": fnum(s.get("perm_P")),
                "fig5a_conc_rate": fnum(s.get("conc_rate")),
                "fig5a_bg_rate": fnum(s.get("bg_rate")),
                "fig5a_selected": sstr(s.get("sel")),
            }
        traits.append({
            "id": tid,
            "short": sstr(sub.get("short")),
            "source": sstr(sub.get("source")),
            "mask": sstr(sub.get("mask")),
            "class": sstr(sub.get("program_class")),
            "flag": sstr(tr.get("flag")) or "ok",
            "desc": cfg["trait_desc"].get(tid, ""),
            "lof_thresh": fnum(tr.get("thresh")),
            "n_top": fnum(tr.get("n_top")),
            "pct_neg": fnum(tr.get("pct_neg")),
            "r_prior": fnum(tr.get("r_prior")),
            "panels": panels,
        })
    short_by_trait = {t["id"]: (t["short"] or t["id"]) for t in traits}
    flag_by_trait = {t["id"]: t["flag"] for t in traits}

    # ---- programs per condition ----
    meta_idx = {(r["program"], r["cond"]): r for _, r in meta.iterrows()}
    programs = {"rest": [], "stim48hr": []}
    for cond in ("rest", "stim48hr"):
        tc = tidy[tidy["cond"] == cond]
        for prog in sorted(tc["program"].unique(), key=prog_num):
            rows = tc[tc["program"] == prog]
            ann = ""
            for a in rows["annotation"]:
                if sstr(a):
                    ann = sstr(a); break
            per_trait = {}
            for _, r in rows.iterrows():
                tid = r["trait"]
                per_trait[short_by_trait.get(tid, tid)] = {
                    "trait_id": tid,
                    "source": sstr(r.get("source")),
                    "flag": flag_by_trait.get(tid, "ok"),
                    "prog_P": fnum(r.get("prog_P")),
                    "prog_FDR": fnum(r.get("prog_FDR")),
                    "prog_dir": fnum(r.get("prog_dir")),
                    "prog_effect": fnum(r.get("prog_effect")),
                    "reg_P": fnum(r.get("reg_P_withShet")),
                    "reg_FDR": fnum(r.get("reg_FDR")),
                    "reg_dir": fnum(r.get("reg_dir")),
                    "reg_beta": fnum(r.get("reg_beta_withShet")),
                    "reg_beta_se": fnum(r.get("reg_betaSE_withShet")),
                }
            mrow = meta_idx.get((prog, cond))
            m = None
            if mrow is not None:
                m = {
                    "pooled": fnum(mrow.get("all5_pooled")),
                    "pooled_se": fnum(mrow.get("all5_pooled_se")),
                    "P": fnum(mrow.get("all5_P")),
                    "FDR": fnum(mrow.get("all5_FDR")),
                    "I2": fnum(mrow.get("all5_I2")),
                    "k": fnum(mrow.get("all5_k")),
                    "n_pos": fnum(mrow.get("all5_n_pos")),
                    "n_neg": fnum(mrow.get("all5_n_neg")),
                    "loo_sign_stable": sstr(mrow.get("loo_sign_stable")) == "True",
                    "cross_source_concordant": sstr(mrow.get("cross_source_concordant")) == "True",
                    "backman_pooled": fnum(mrow.get("Backman_pooled")),
                    "genebass_pooled": fnum(mrow.get("GeneBass_pooled")),
                }
            # a compact "best evidence" number for default sorting in the UI
            best_prog_P = min([v["prog_P"] for v in per_trait.values() if v["prog_P"] is not None],
                              default=None)
            programs[cond].append({
                "program": prog,
                "num": prog_num(prog),
                "annotation": ann,
                "top_genes": topg.get(cond, {}).get(prog_num(prog), []),
                "meta": m,
                "best_prog_P": best_prog_P,
                "meta_FDR": (m or {}).get("FDR"),
                "traits": per_trait,
            })

    figures = copy_figures(did, cfg, [t["id"] for t in traits])
    narrative = load_narrative(rr)

    # ---- lookups for the interactive network + regulator modules ----
    ann_by = {c: {p["program"]: p["annotation"] for p in programs[c]} for c in programs}
    meta_sig = {c: {p["program"] for p in programs[c]
                    if (p.get("meta") or {}).get("FDR") is not None and p["meta"]["FDR"] < 0.05}
                for c in programs}
    progMeta_by = {}
    for c in programs:
        for p in programs[c]:
            for short, tv in p["traits"].items():
                progMeta_by[(c, p["program"], short)] = {
                    "P": tv.get("prog_P"), "dir": tv.get("prog_dir"), "effect": tv.get("prog_effect")}
    flag_by_short = {t["short"]: t["flag"] for t in traits}

    networks = build_networks(rr, cfg, ann_by, meta_sig, progMeta_by, short_by_trait)
    regulators = build_regulators(rr, [t["id"] for t in traits], ann_by, short_by_trait, flag_by_short)

    return {
        "id": did,
        "name": cfg["name"],
        "cell_type": cfg["cell_type"],
        "conditions": ["rest", "stim48hr"],
        "n_programs": 60,
        "traits": traits,
        "programs": programs,
        "figures": figures,
        "networks": networks,
        "regulators": regulators,
        "narrative": narrative,
    }

# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def _caption(kind, trait_short, cond):
    if kind == "Fig4C":
        return f"Fig 4C — program vs regulator burden ({trait_short}, {cond})"
    if kind == "Fig5a":
        return f"Fig 5a — regulators→programs→trait map ({trait_short}, {cond})"
    return kind

def copy_figures(did, cfg, trait_ids):
    rr = cfg["results_root"]
    out_dir = os.path.join(HERE, "static", "figures", did)
    os.makedirs(out_dir, exist_ok=True)
    figs = []

    def take(src, meta):
        if not os.path.exists(src):
            return
        base = os.path.basename(src)
        dst = os.path.join(out_dir, base)
        shutil.copyfile(src, dst)
        meta["file"] = f"figures/{did}/{base}"
        figs.append(meta)

    short_map = {}
    for tid in trait_ids:
        short_map[tid] = tid.replace("_conv", "").replace("Backman_2021_", "BK_").replace("Genebass_", "GB_")

    # per-trait Fig 4C + Fig 5a map, both conditions
    for tid in trait_ids:
        fdir = os.path.join(rr, tid, "figures")
        for cond in ("rest", "stim48hr"):
            pref = f"{cond}__"
            for f4 in glob.glob(os.path.join(fdir, f"{pref}Fig4C_*.png")):
                take(f4, {"scope": "trait", "trait": tid, "condition": cond, "kind": "Fig4C",
                          "caption": _caption("Fig4C", short_map[tid], cond)})
            for f5 in glob.glob(os.path.join(fdir, f"{pref}*Fig5a_map.png")):
                take(f5, {"scope": "trait", "trait": tid, "condition": cond, "kind": "Fig5a",
                          "caption": _caption("Fig5a", short_map[tid], cond)})

    # global condition-comparison figures
    cc = {
        "1_program_heatmap.png": "Program × trait burden heatmap (per condition)",
        "1b_condition_change_heatmap.png": "Rest→Stim change (Δβ), significance + direction",
        "2_rest_vs_stim_scatter.png": "Rest vs Stim48hr regulator β (matched programs)",
        "3_forest_priority.png": "Random-effects meta-analysis forest (priority programs)",
        "4_sign_overlap.png": "Sign-concordance overlap vs permutation null",
        "5_p58_qc.png": "P58 olfactory/noise QC exclusion",
    }
    for fn, cap in cc.items():
        take(os.path.join(rr, "condition_comparison/figures", fn),
             {"scope": "global", "trait": None, "condition": None, "kind": "condition_comparison",
              "caption": cap})

    # combined annotated panels
    for src in sorted(glob.glob(os.path.join(rr, "combined_figures_annotated", "*.png"))):
        base = os.path.basename(src)
        cond = "stim48hr" if "stim48" in base else "rest"
        grp = "Group1 RA" if "Group1" in base else ("Group2 RF" if "Group2" in base else "")
        take(src, {"scope": "global", "trait": None, "condition": cond, "kind": "combined",
                   "caption": f"Combined Fig 5a (annotated) — {grp}, {cond}"})
    return figs

def load_narrative(rr):
    out = {}
    for key, path in (("conclusions", "condition_comparison/CONCLUSIONS.md"),
                      ("scope", "condition_comparison/DATA_AUDIT_AND_SCOPE.md"),
                      ("overview", "README.md")):
        fp = os.path.join(rr, path)
        if os.path.exists(fp):
            with open(fp) as fh:
                out[key] = fh.read()
    return out

# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--disease", default="RA", help="disease id present in DISEASES")
    ap.add_argument("--all", action="store_true", help="build every disease in DISEASES")
    args = ap.parse_args()

    ids = list(DISEASES) if args.all else [args.disease]
    data = {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "diseases": {}}
    for did in ids:
        if did not in DISEASES:
            sys.exit(f"unknown disease '{did}'; known: {list(DISEASES)}")
        print(f"[build] {did} ...", file=sys.stderr)
        data["diseases"][did] = build_disease(did, DISEASES[did])

    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    out = os.path.join(HERE, "data", "demo_data.json")
    with open(out, "w") as fh:
        json.dump(data, fh, allow_nan=False, separators=(",", ":"))
    # report
    for did, d in data["diseases"].items():
        nfig = len(d["figures"])
        print(f"[ok] {did}: {len(d['traits'])} traits, "
              f"{len(d['programs']['rest'])}/{len(d['programs']['stim48hr'])} programs (rest/stim), "
              f"{nfig} figures", file=sys.stderr)
    print(f"[wrote] {out} ({os.path.getsize(out)//1024} KB)", file=sys.stderr)

if __name__ == "__main__":
    main()
