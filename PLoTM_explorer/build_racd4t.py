#!/usr/bin/env python
"""build_racd4t.py — flatten the RA_CD4T_external analysis (patient CD4+ T cNMF, K=30, projected onto
Perturb-seq in 2 states, burden-tested against 3 GeneBass traits) into demo_data.json, in the SAME
schema build_demo_data.py emits, so the explorer + network app surface it as the default dataset.

    PYTHONNOUSERSITE=1 /mnt/scratch/ZY2/.envs/perturbseq-repro/bin/python build_racd4t.py

Source: tcell_perturbseq/cNMF_RA_analysis/RA_CD4T_external/
Unlike the Perturb-seq RA run, the cNMF here is a SINGLE patient-derived decomposition (so a program has
the same identity in both states); only the regulatory betas and burden differ by state.
"""
import glob
import json
import math
import os
import shutil
import sys
from datetime import datetime, timezone

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("PROGLOF_ROOT", os.path.abspath(os.path.join(HERE, "..")))
SRC = os.path.join(ROOT, "tcell_perturbseq/cNMF_RA_analysis/RA_CD4T_external")
DATA = os.path.join(ROOT, "Hackthron-Claude/Paper/reproduction/data")   # gamma / gencode / shet
REGDIR = os.path.join(SRC, "stageBC/run_root/data/Perturbseq/cNMF_regulation")
DID = "RA_CD4T"
FIGDIR = os.path.join(HERE, "static", "figures", DID)

TRAITS = [  # (full id, short, class, description)
    ("Genebass_RA_custom", "GB_RA_custom", "binary", "RA, custom case definition"),
    ("Genebass_RA_M06", "GB_RA_M06", "binary", "RA, phecode M06"),
    ("Genebass_RF", "GB_RF", "continuous", "rheumatoid factor (endophenotype)"),
]
SHORT = {t[0]: t[1] for t in TRAITS}
NET_TRAITS = ["Genebass_RA_custom", "Genebass_RA_M06"]  # the traits with Fig-5 maps
STATES = [("rest", "rest"), ("stim48", "stim48hr")]      # (source state, canonical condition)
CANON = {"rest": "rest", "stim48": "stim48hr"}


def fnum(x):
    try:
        v = float(x)
        return None if (math.isnan(v) or math.isinf(v)) else v
    except (TypeError, ValueError):
        return None


def bh_fdr(pmap):
    """Benjamini-Hochberg over {key: p}; returns {key: q}."""
    items = [(k, p) for k, p in pmap.items() if p is not None]
    n = len(items)
    if not n:
        return {}
    items.sort(key=lambda kv: kv[1])
    q = {}
    prev = 1.0
    for i in range(n - 1, -1, -1):
        k, p = items[i]
        prev = min(prev, p * n / (i + 1))
        q[k] = prev
    return q


def meta_dl(betas, ses):
    """DerSimonian-Laird random-effects meta of beta ± se (drop missing)."""
    pairs = [(b, s) for b, s in zip(betas, ses) if b is not None and s is not None and s > 0]
    k = len(pairs)
    if k == 0:
        return None
    w = [1 / s ** 2 for _, s in pairs]
    fe = sum(wi * b for wi, (b, _) in zip(w, pairs)) / sum(w)
    Q = sum(wi * (b - fe) ** 2 for wi, (b, _) in zip(w, pairs))
    df = k - 1
    C = sum(w) - sum(wi ** 2 for wi in w) / sum(w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    wr = [1 / (s ** 2 + tau2) for _, s in pairs]
    pooled = sum(wri * b for wri, (b, _) in zip(wr, pairs)) / sum(wr)
    se = math.sqrt(1 / sum(wr))
    z = pooled / se if se > 0 else 0.0
    P = math.erfc(abs(z) / math.sqrt(2))
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    npos = sum(1 for b, _ in pairs if b > 0)
    # leave-one-out sign stability
    stable = True
    if k >= 2:
        full_sign = 1 if pooled >= 0 else -1
        for j in range(k):
            sub = pairs[:j] + pairs[j + 1:]
            ww = [1 / s ** 2 for _, s in sub]
            p = sum(wi * b for wi, (b, _) in zip(ww, sub)) / sum(ww)
            if (1 if p >= 0 else -1) != full_sign:
                stable = False
                break
    return {"pooled": pooled, "pooled_se": se, "P": P, "I2": I2, "k": k,
            "n_pos": npos, "n_neg": k - npos, "loo_sign_stable": stable,
            "cross_source_concordant": None, "backman_pooled": None, "genebass_pooled": pooled}


# --------------------------------------------------------------------------- readers
def read_annotation():
    df = pd.read_csv(os.path.join(SRC, "annotation/program_annotation_named.tsv"), sep="\t", dtype=str)
    out = {}
    for _, r in df.iterrows():
        out[f"P{int(r['Program'])}"] = {"name": r.get("name", ""), "kind": r.get("kind", "")}
    return out


def read_top_genes(topn=12):
    p = os.path.join(SRC, "cNMF/RA_CD4T/test1/test1.gene_spectra_score.k_30.dt_0_4.txt")
    df = pd.read_csv(p, sep="\t", index_col=0)
    out = {}
    for prog, row in df.iterrows():
        out[f"P{int(prog)}"] = list(row.sort_values(ascending=False).head(topn).index)
    return out


def read_burden():
    """{(pid, state, trait): {prog_*, reg_*}}"""
    df = pd.read_csv(os.path.join(SRC, "stageBC/RA_program_x_2states_x_3traits_burden.tsv"),
                     sep="\t", dtype=str)
    out = {}
    for _, r in df.iterrows():
        pid = f"P{int(r['Program'])}"
        meang = fnum(r.get("MEANgamma_top100"))
        rand = fnum(r.get("shet_adjusted_random_mean"))
        eff = (meang - rand) if (meang is not None and rand is not None) else None
        beta = fnum(r.get("beta_withShet"))
        out[(pid, r["state"], r["trait"])] = {
            "prog_P": fnum(r.get("MEANgamma_top100_shet_adjusted_P")),
            "prog_effect": eff, "prog_dir": (1 if (eff or 0) >= 0 else -1),
            "reg_P": fnum(r.get("P_withShet")), "reg_beta": beta,
            "reg_beta_se": fnum(r.get("betaSE_withShet")),
            "reg_dir": (1 if (beta or 0) >= 0 else -1),
            "estimable": r.get("estimable") == "True",
        }
    return out


def read_fig5_summary():
    p = os.path.join(SRC, "fig5/tables/fig5_summary_RAprog.tsv")
    df = pd.read_csv(p, sep="\t", dtype=str)
    out = {}
    for _, r in df.iterrows():
        out[(r["DS"], r["TRAIT"])] = {"perm_P": fnum(r.get("perm_P")),
                                      "conc_rate": fnum(r.get("conc_rate")),
                                      "bg_rate": fnum(r.get("bg_conc_rate"))}
    return out


def _map_paths(state, trait):
    base = os.path.join(SRC, f"fig5/out/RAprog_GWCD4i_{state}/P200")
    tag = f"{trait}_program5_regulator3_LOF0.03"
    return (os.path.join(base, f"{tag}_map_nodes.tsv"),
            os.path.join(base, f"{tag}_map_edges.tsv"),
            os.path.join(base, f"SelectedPrograms_{trait}_program5_regulator3_LOF0.03.txt"))


def _pid(x):
    x = str(x)
    return x if x.upper().startswith("P") else f"P{x}"


# --------------------------------------------------------------------------- builders
def build_programs(ann, topg, burden, meta_sig_ref):
    programs = {"rest": [], "stim48hr": []}
    meta_all = {"rest": {}, "stim48hr": {}}
    # FDR per (trait, condition)
    fdr_prog, fdr_reg = {}, {}
    for state, cond in STATES:
        for tid, *_ in TRAITS:
            pm = {f"P{i}": burden.get((f"P{i}", state, tid), {}).get("prog_P") for i in range(1, 31)}
            rm = {f"P{i}": burden.get((f"P{i}", state, tid), {}).get("reg_P") for i in range(1, 31)}
            fdr_prog[(cond, tid)] = bh_fdr(pm)
            fdr_reg[(cond, tid)] = bh_fdr(rm)
    for state, cond in STATES:
        for i in range(1, 31):
            pid = f"P{i}"
            a = ann.get(pid, {})
            betas = [burden.get((pid, state, tid), {}).get("reg_beta") for tid, *_ in TRAITS]
            ses = [burden.get((pid, state, tid), {}).get("reg_beta_se") for tid, *_ in TRAITS]
            meta = meta_dl(betas, ses)
            per_trait = {}
            for tid, short, *_ in TRAITS:
                b = burden.get((pid, state, tid), {})
                per_trait[short] = {
                    "trait_id": tid, "source": "GeneBass",
                    "flag": "ok" if b.get("estimable", True) else "not-estimable",
                    "prog_P": b.get("prog_P"), "prog_FDR": fdr_prog[(cond, tid)].get(pid),
                    "prog_dir": b.get("prog_dir"), "prog_effect": b.get("prog_effect"),
                    "reg_P": b.get("reg_P"), "reg_FDR": fdr_reg[(cond, tid)].get(pid),
                    "reg_dir": b.get("reg_dir"), "reg_beta": b.get("reg_beta"),
                    "reg_beta_se": b.get("reg_beta_se"),
                }
            best = min([v["prog_P"] for v in per_trait.values() if v["prog_P"] is not None], default=None)
            programs[cond].append({
                "program": pid, "num": i,
                "annotation": a.get("name", "") if a.get("kind") == "biological" else a.get("name", ""),
                "kind": a.get("kind", ""), "top_genes": topg.get(pid, []),
                "meta": meta, "best_prog_P": best,
                "meta_FDR": (meta or {}).get("FDR"), "traits": per_trait,
            })
        # meta FDR across programs per condition
        mp = {p["program"]: (p["meta"] or {}).get("P") for p in programs[cond]}
        mq = bh_fdr(mp)
        for p in programs[cond]:
            if p["meta"] is not None:
                p["meta"]["FDR"] = mq.get(p["program"])
                p["meta_FDR"] = mq.get(p["program"])
            if (p["meta"] or {}).get("FDR") is not None and p["meta"]["FDR"] < 0.05:
                meta_sig_ref.setdefault(cond, set()).add(p["program"])
    return programs


def build_networks(ann, burden, meta_sig):
    label = "RA programs (custom + M06)"
    out = {"rest": {}, "stim48hr": {}}
    for state, cond in STATES:
        gene_sign, gene_traits, prog_traits, prog_role, edge_agg = {}, {}, {}, {}, {}
        trait_ids = []
        for tid in NET_TRAITS:
            short = SHORT[tid]
            trait_ids.append(short)
            nfile, efile, sfile = _map_paths(state, tid)
            if not os.path.exists(nfile):
                continue
            nodes = pd.read_csv(nfile, sep="\t", dtype=str)
            edges = pd.read_csv(efile, sep="\t", dtype=str) if os.path.getsize(efile) else pd.DataFrame(columns=["GENE", "PROGRAM", "beta", "padj", "type"])
            psel, rsel = [], []
            if os.path.exists(sfile):
                sdf = pd.read_csv(sfile, sep="\t", dtype=str)
                if len(sdf):
                    psel = [_pid(x) for x in str(sdf.iloc[0].get("Program_selected", "")).split(",") if x]
                    rsel = [_pid(x) for x in str(sdf.iloc[0].get("Regulator_selected", "")).split(",") if x]
            for _, r in nodes.iterrows():
                if r["kind"] == "gene":
                    gene_sign.setdefault(r["id"], []).append(int(fnum(r.get("sign")) or 0))
                    gene_traits.setdefault(r["id"], set()).add(short)
                elif r["kind"] == "program":
                    prog_traits.setdefault(_pid(r["id"]), set()).add(short)
            for p in psel:
                prog_role.setdefault(p, set()).add("program"); prog_traits.setdefault(p, set()).add(short)
            for p in rsel:
                prog_role.setdefault(p, set()).add("regulator"); prog_traits.setdefault(p, set()).add(short)
            for _, r in edges.iterrows():
                key = (r["GENE"], _pid(r["PROGRAM"]))
                e = edge_agg.setdefault(key, {"type": r.get("type"), "betas": [], "traits": set()})
                b = fnum(r.get("beta"))
                if b is not None:
                    e["betas"].append(b)
                e["traits"].add(short)
                if r.get("type") == "regulates":
                    e["type"] = "regulates"
                gene_traits.setdefault(r["GENE"], set()).add(short)
        gene_is_reg = {g for (g, _), e in edge_agg.items() if e["type"] == "regulates"}
        nodes_out, edges_out = [], []
        for gene, traits in gene_traits.items():
            signs = [s for s in gene_sign.get(gene, []) if s]
            sign = max(set(signs), key=signs.count) if signs else 0
            nodes_out.append({"id": gene, "type": "gene",
                              "role": "regulator" if gene in gene_is_reg else "member",
                              "sign": sign, "n_traits": len(traits), "traits": sorted(traits)})
        for prog, traits in prog_traits.items():
            roles = prog_role.get(prog, set())
            role = "both" if len(roles) == 2 else (next(iter(roles)) if roles else "program")
            nm = ann.get(prog, {}).get("name", "")
            nodes_out.append({"id": prog, "type": "program", "annotation": nm,
                              "annotated": ann.get(prog, {}).get("kind") == "biological", "role": role,
                              "n_traits": len(traits), "traits": sorted(traits),
                              "sig": prog in meta_sig.get(cond, set())})
        for tn in sorted(set(trait_ids)):
            nodes_out.append({"id": tn, "type": "trait"})
        for (gene, prog), e in edge_agg.items():
            beta = max(e["betas"], key=abs) if e["betas"] else None
            edges_out.append({"source": gene, "target": prog, "cls": e["type"], "beta": beta,
                              "n_traits": len(e["traits"]), "traits": sorted(e["traits"])})
        for prog, traits in prog_traits.items():
            for short in traits:
                tid = next(t[0] for t in TRAITS if t[1] == short)
                b = burden.get((prog, state, tid), {})
                w = b.get("prog_P")
                edges_out.append({"source": prog, "target": short, "cls": "selects",
                                  "weight": (-math.log10(w) if (w and w > 0) else 1.0),
                                  "dir": b.get("prog_dir", -1), "effect": b.get("prog_effect")})
        out[cond]["Group_RA"] = {"label": label, "nodes": nodes_out, "edges": edges_out,
                                 "n_programs": sum(1 for n in nodes_out if n["type"] == "program"),
                                 "n_genes": sum(1 for n in nodes_out if n["type"] == "gene")}
    return out


def _gencode():
    df = pd.read_csv(os.path.join(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID"), sep="\t",
                     header=None, names=["ensg", "sym"], dtype=str)
    df = df.drop_duplicates("ensg")
    return dict(zip(df.ensg, df.sym))


def _read_gamma(trait_file, ens2sym):
    """{symbol: gamma}, {symbol: ensg} from the GeneBayes posterior."""
    p = os.path.join(DATA, "LoF/GeneBayes_posterior", trait_file)
    df = pd.read_csv(p, sep="\t", dtype=str)
    g, e = {}, {}
    for _, r in df.iterrows():
        sym = ens2sym.get(r["ensg"])
        v = fnum(r.get("post_mean"))
        if sym and v is not None and sym not in g:
            g[sym] = v
            e[sym] = r["ensg"]
    return g, e


def _read_shet():
    df = pd.read_csv(os.path.join(DATA, "shet_10bins.txt"), sep="\t", dtype=str)
    return {r["ensg"]: fnum(r["shet"]) for _, r in df.iterrows()}


def _read_perturb(state, pid):
    """{gene: (beta, p)} for a regulator program in a state."""
    i = int(str(pid).lstrip("Pp"))
    f = os.path.join(REGDIR, f"RAprog_GWCD4i_{state}", f"K30_program{i}_perturb_effects.txt")
    if not os.path.exists(f):
        return {}
    df = pd.read_csv(f, sep="\t", dtype=str)
    return {r["GENE"]: (fnum(r["lm_es"]), fnum(r["lm_p"])) for _, r in df.iterrows()}


def _bh(pvals):
    import numpy as np
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for k in range(n - 1, -1, -1):
        idx = order[k]
        prev = min(prev, p[idx] * n / (k + 1))
        q[idx] = prev
    return q


def _joint_wp(gamma, ens_of, shet, betas_by_prog, reg):
    """Reproduce fig5a_regulators_RAprog.R: scaled OLS  gamma ~ sum_P beta_P + shet ; return
    {P: (wP, wP_pvalue)} for the regulator programs."""
    import numpy as np
    genes = set(gamma) & set(ens_of)
    for p in reg:
        genes &= set(betas_by_prog[p])
    genes = [g for g in genes if ens_of[g] in shet and shet[ens_of[g]] is not None]
    if len(genes) < len(reg) + 3:
        return {p: (None, None) for p in reg}
    y = np.array([gamma[g] for g in genes], float)
    cols = [np.array([betas_by_prog[p][g][0] for g in genes], float) for p in reg]
    sv = np.array([shet[ens_of[g]] for g in genes], float)
    mats = [y] + cols + [sv]
    scaled = []
    for v in mats:
        v = v.copy()
        fin = np.isfinite(v)
        if not fin.all():
            v[~fin] = np.nanmax(v[fin]) if fin.any() else 0.0
        sd = v.std(ddof=1)
        scaled.append((v - v.mean()) / sd if sd > 0 else v * 0.0)
    ys = scaled[0]
    X = np.column_stack([np.ones(len(genes))] + scaled[1:])  # intercept + betas + shet
    coef, *_ = np.linalg.lstsq(X, ys, rcond=None)
    resid = ys - X @ coef
    dof = len(genes) - X.shape[1]
    sigma2 = float(resid @ resid) / dof if dof > 0 else np.nan
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(sigma2 * XtX_inv))
    try:
        from scipy import stats
        pv = 2 * stats.t.sf(np.abs(coef / se), dof)
    except Exception:  # noqa: BLE001 — normal approximation
        pv = np.array([math.erfc(abs(z) / math.sqrt(2)) for z in (coef / se)])
    return {p: (float(coef[1 + j]), float(pv[1 + j])) for j, p in enumerate(reg)}


def build_reg_networks(ann):
    """Paper's Fig-5a 'regulator model', COMBINED across the two RA traits per condition: the union of
    each trait's RN regsubsets-selected regulator programs, each with its top-10 trans-regulators
    (BH<0.05, by |beta|); gene→program coloured by sign(beta); each program links to the trait(s) that
    selected it with THAT trait's joint-regression weight w_P; two trait nodes on the right."""
    ens2sym = _gencode()
    shet = _read_shet()
    out = {"rest": {}, "stim48hr": {}}
    for state, cond in STATES:
        prog_nodes = {}          # pid -> {annotation, annotated, traits:set, sig}
        greg = {}                # (gene, pid) -> beta   (trait-independent regulation; deduped)
        gene_gamma = {}          # gene -> [gamma per trait it regulates for]
        gene_traits = {}         # gene -> set of trait shorts
        pt_edges = []            # program -> trait edges (per trait, with that trait's w_P)
        trait_shorts = []
        for tid in NET_TRAITS:
            short = SHORT[tid]
            _, _, sfile = _map_paths(state, tid)
            if not os.path.exists(sfile):
                continue
            trait_shorts.append(short)
            sdf = pd.read_csv(sfile, sep="\t", dtype=str)
            reg = [_pid(x) for x in str(sdf.iloc[0].get("Regulator_selected", "")).split(",") if x]
            betas_by_prog = {p: _read_perturb(state, p) for p in reg}
            gamma, ens_of = _read_gamma(f"{tid}.per_gene_estimates.tsv", ens2sym)
            wp = _joint_wp(gamma, ens_of, shet, betas_by_prog, reg)
            for p in reg:
                be = betas_by_prog[p]
                genes = list(be.keys())
                pv = [be[g][1] if be[g][1] is not None else 1.0 for g in genes]
                q = _bh(pv)
                top = [(g, be[g][0]) for g, qq in zip(genes, q) if qq < 0.05 and be[g][0] is not None]
                top.sort(key=lambda gb: -abs(gb[1]))
                for g, b in top[:10]:
                    greg[(g, p)] = b  # same regulation regardless of trait -> dedupe
                    gene_gamma.setdefault(g, []).append(gamma.get(g))
                    gene_traits.setdefault(g, set()).add(short)
                pn = prog_nodes.setdefault(p, {"annotation": ann.get(p, {}).get("name", ""),
                                               "annotated": ann.get(p, {}).get("kind") == "biological",
                                               "traits": set(), "sig": False})
                pn["traits"].add(short)
                w, wpp = wp.get(p, (None, None))
                if wpp is not None and wpp < 0.05:
                    pn["sig"] = True
                if w is not None:
                    pt_edges.append({"source": p, "target": short, "cls": "selects", "wP": w, "wP_P": wpp,
                                     "dir": (1 if w > 0 else -1),
                                     "weight": (-math.log10(wpp) if (wpp and wpp > 0) else 1.0),
                                     "label": f"w={w:+.3f}\nP={wpp:.2g}"})
        nodes, edges = [], []
        for (g, p), b in greg.items():
            edges.append({"source": g, "target": p, "cls": "regulates", "beta": b,
                          "dir": (1 if b > 0 else -1), "n_traits": len(gene_traits.get(g, [])),
                          "traits": sorted(gene_traits.get(g, []))})
        edges += pt_edges
        for g in sorted(gene_gamma):
            vals = [x for x in gene_gamma[g] if x is not None]
            gv = (sum(vals) / len(vals)) if vals else None
            gsign = 0 if gv is None else (1 if gv > 0 else (-1 if gv < 0 else 0))
            gclass = "small" if (gv is None or abs(gv) <= 0.03) else ("pos" if gv > 0 else "neg")
            tr = sorted(gene_traits.get(g, []))
            nodes.append({"id": g, "type": "gene", "role": "regulator", "gamma": gv,
                          "sign": gsign, "gclass": gclass, "n_traits": len(tr), "traits": tr})
        for p, pn in prog_nodes.items():
            nodes.append({"id": p, "type": "program", "annotation": pn["annotation"],
                          "annotated": pn["annotated"], "role": "regulator",
                          "n_traits": len(pn["traits"]), "traits": sorted(pn["traits"]), "sig": pn["sig"]})
        for short in dict.fromkeys(trait_shorts):
            nodes.append({"id": short, "type": "trait"})
        out[cond]["combined"] = {"label": "RA_custom + RA_M06 — regulator model", "trait": None,
                                 "nodes": nodes, "edges": edges,
                                 "n_programs": len(prog_nodes), "n_genes": len(gene_gamma)}
    return out


def build_regulators(ann):
    reg, gene_sign = {}, {}
    for state, cond in STATES:
        for tid in NET_TRAITS:
            short = SHORT[tid]
            nfile, efile, _ = _map_paths(state, tid)
            if not os.path.exists(efile):
                continue
            if os.path.exists(nfile):
                nodes = pd.read_csv(nfile, sep="\t", dtype=str)
                for _, r in nodes.iterrows():
                    if r["kind"] == "gene":
                        gene_sign.setdefault(r["id"], []).append(int(fnum(r.get("sign")) or 0))
            if not os.path.getsize(efile):
                continue
            edges = pd.read_csv(efile, sep="\t", dtype=str)
            for _, r in edges.iterrows():
                if r.get("type") != "regulates":
                    continue
                gene, prog, b = r["GENE"], _pid(r["PROGRAM"]), fnum(r.get("beta"))
                g = reg.setdefault(gene, {"targets": [], "traits": set(), "conditions": set(), "programs": set()})
                g["targets"].append({"program": prog, "annotation": ann.get(prog, {}).get("name", ""),
                                     "beta": b, "dir": (1 if (b or 0) > 0 else -1),
                                     "trait": short, "condition": cond, "flag": "ok"})
                g["traits"].add(short); g["conditions"].add(cond); g["programs"].add(prog)
    out = []
    for gene, g in reg.items():
        signs = [s for s in gene_sign.get(gene, []) if s]
        sign = max(set(signs), key=signs.count) if signs else 0
        ct = sorted(g["traits"])
        out.append({"gene": gene, "sign": sign, "n_programs": len(g["programs"]),
                    "n_targets": len(g["targets"]), "n_traits": len(g["traits"]),
                    "traits": ct, "conditions": sorted(g["conditions"]),
                    "n_clean_traits": len(ct), "clean_traits": ct, "shared": len(ct) >= 2,
                    "targets": sorted(g["targets"], key=lambda t: -abs(t["beta"] or 0))})
    out.sort(key=lambda r: (-r["n_clean_traits"], -r["n_programs"], r["gene"]))
    return out


def copy_figures(fig5_summary):
    os.makedirs(FIGDIR, exist_ok=True)
    figs = []

    def take(src, meta):
        if not os.path.exists(src):
            return
        base = os.path.basename(src)
        shutil.copyfile(src, os.path.join(FIGDIR, base))
        meta["file"] = f"figures/{DID}/{base}"
        figs.append(meta)

    for state, cond in STATES:
        for tid, short, *_ in TRAITS:
            take(os.path.join(SRC, f"stageBC/figures_{state}/Fig4C_RA_programs_GWCD4i_{state}_{tid}.png"),
                 {"scope": "trait", "trait": tid, "condition": cond, "kind": "Fig4C",
                  "caption": f"Fig 4C — program vs regulator burden ({short}, {cond})"})
        for tid in NET_TRAITS:
            short = SHORT[tid]
            take(os.path.join(SRC, f"fig5/figures/{tid}_program5_regulator3_LOF0.03_RAprog_GWCD4i_{state}_Fig5a_map.png"),
                 {"scope": "trait", "trait": tid, "condition": cond, "kind": "Fig5a",
                  "caption": f"Fig 5a — regulators→programs→trait ({short}, {cond})"})
        take(os.path.join(SRC, f"fig5/figures/Fig5a_twoTrait_RAprog_GWCD4i_{state}.png"),
             {"scope": "global", "trait": None, "condition": cond, "kind": "combined",
              "caption": f"Combined Fig 5a — RA custom + M06 ({cond})"})
    return figs


def load_narrative():
    out = {}
    for key, path in (("overview", "README.md"), ("scope", "stageBC/README.md"),
                      ("conclusions", "fig5/BIOLOGY_NOTES.md")):
        fp = os.path.join(SRC, path)
        if os.path.exists(fp):
            out[key] = open(fp).read()
    return out


def main():
    ann = read_annotation()
    topg = read_top_genes()
    burden = read_burden()
    fig5 = read_fig5_summary()

    meta_sig = {}
    programs = build_programs(ann, topg, burden, meta_sig)

    # traits[] with per-condition panels
    traits = []
    for tid, short, cls, desc in TRAITS:
        panels = {}
        for state, cond in STATES:
            s = fig5.get((state, tid))
            # Fig 4C Bonferroni hit = program with smallest prog_P < 0.05/30
            hits = sorted([(burden.get((f"P{i}", state, tid), {}).get("prog_P"), f"P{i}") for i in range(1, 31)
                           if burden.get((f"P{i}", state, tid), {}).get("prog_P") is not None])
            hit = hits[0][1] if hits and hits[0][0] is not None and hits[0][0] < 0.05 / 30 else "-"
            panels[cond] = {"fig4c_prog_hit": hit, "fig4c_reg_hit": "-",
                            "fig5a_perm_P": (s or {}).get("perm_P"),
                            "fig5a_conc_rate": (s or {}).get("conc_rate"),
                            "fig5a_bg_rate": (s or {}).get("bg_rate"), "fig5a_selected": ""}
        traits.append({"id": tid, "short": short, "source": "GeneBass", "mask": "M1_pLoF",
                       "class": cls, "flag": "ok", "desc": desc, "lof_thresh": 0.03,
                       "n_top": 200, "pct_neg": None, "r_prior": None, "panels": panels})

    disease = {
        "id": DID, "name": "Rheumatoid arthritis — patient CD4⁺ T programs",
        "cell_type": "Patient CD4⁺ T cells (external RA scRNA-seq) · cNMF K=30 · projected onto CRISPRi Perturb-seq",
        "conditions": ["rest", "stim48hr"], "n_programs": 30,
        "traits": traits, "programs": programs,
        "figures": copy_figures(fig5),
        "networks": build_reg_networks(ann),
        "regulators": build_regulators(ann),
        "narrative": load_narrative(),
    }
    data = {"generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "diseases": {DID: disease}}
    os.makedirs(os.path.join(HERE, "data"), exist_ok=True)
    out = os.path.join(HERE, "data", "demo_data.json")
    with open(out, "w") as fh:
        json.dump(data, fh, allow_nan=False, separators=(",", ":"))
    nfig = len(disease["figures"])
    print(f"[ok] {DID}: {len(traits)} traits, 30 programs x 2 conds, {nfig} figures, "
          f"{len(disease['regulators'])} regulators", file=sys.stderr)
    print(f"[wrote] {out} ({os.path.getsize(out)//1024} KB)", file=sys.stderr)


if __name__ == "__main__":
    main()
