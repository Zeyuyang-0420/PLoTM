"""proglof_data.py — read-only accessor over data/demo_data.json.

Dependency-free (json + os only) so it can be imported by the FastAPI chat backend AND, unchanged, by
the Claude-Science `proglof-mcp` connector's read tools — one results-access contract for both
deliverables. Every function returns plain JSON-serialisable dicts/lists.

Vocabulary:
  disease   e.g. "RA"
  condition "rest" | "stim48hr"
  program   "P1".."P60"     (cNMF is independent per condition, so rest-P2 != stim-P2)
  trait     short id, e.g. "GB_RF", "BK_curated"   (8 RA LoF phenotypes)
Scientific grounding rules the callers must honour live in system_brief() below.
"""
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA_PATH = os.environ.get("PROGLOF_DEMO_DATA", os.path.join(_HERE, "data", "demo_data.json"))
_CACHE = {}


def load(path=None):
    p = path or _DATA_PATH
    if p not in _CACHE:
        with open(p) as fh:
            _CACHE[p] = json.load(fh)
    return _CACHE[p]


def _disease(did):
    d = load()["diseases"]
    if did not in d:
        raise KeyError(f"unknown disease '{did}'; have {list(d)}")
    return d[did]


# --------------------------------------------------------------------------- catalogue
def list_diseases():
    out = []
    for did, d in load()["diseases"].items():
        out.append({"id": did, "name": d["name"], "cell_type": d["cell_type"],
                    "conditions": d["conditions"], "n_traits": len(d["traits"])})
    return out


def list_traits(did):
    return _disease(did)["traits"]


# --------------------------------------------------------------------------- programs
def _prog_index(did, condition):
    return {p["program"]: p for p in _disease(did)["programs"][condition]}


def query_programs(did, condition, trait=None, annotated_only=False, significant=None,
                   sort_by="meta_FDR", limit=25):
    """Return compact program summaries for a condition.

    trait          restrict the per-trait stat shown to this trait short-id (else meta is used)
    annotated_only keep only programs carrying a curated annotation
    significant    None | 'meta' (meta FDR<0.05) | 'program' (any trait program-burden P<0.05)
    sort_by        'meta_FDR' | 'meta_pooled' | 'program' (best per-trait program-burden P) | 'num'
    """
    progs = _disease(did)["programs"][condition]
    rows = []
    for p in progs:
        m = p.get("meta") or {}
        tr = (p["traits"].get(trait) if trait else None)
        row = {
            "program": p["program"], "num": p["num"], "annotation": p["annotation"],
            "top_genes": p["top_genes"][:8],
            "meta_pooled": m.get("pooled"), "meta_FDR": m.get("FDR"), "meta_I2": m.get("I2"),
            "cross_source_concordant": m.get("cross_source_concordant"),
            "loo_sign_stable": m.get("loo_sign_stable"),
            "best_program_P": p.get("best_prog_P"),
        }
        if tr:
            row["trait"] = trait
            row["trait_program_P"] = tr.get("prog_P")
            row["trait_program_FDR"] = tr.get("prog_FDR")
            row["trait_regulator_P"] = tr.get("reg_P")
            row["trait_regulator_beta"] = tr.get("reg_beta")
            row["trait_flag"] = tr.get("flag")
        rows.append((p, row))

    def keep(p, row):
        if annotated_only and not p["annotation"]:
            return False
        if significant == "meta":
            f = (p.get("meta") or {}).get("FDR")
            return f is not None and f < 0.05
        if significant == "program":
            return any((v.get("prog_P") is not None and v["prog_P"] < 0.05)
                       for v in p["traits"].values())
        return True

    rows = [(p, r) for (p, r) in rows if keep(p, r)]

    big = float("inf")
    keyf = {
        "meta_FDR": lambda pr: (pr[1]["meta_FDR"] if pr[1]["meta_FDR"] is not None else big),
        "meta_pooled": lambda pr: -(abs(pr[1]["meta_pooled"]) if pr[1]["meta_pooled"] is not None else -big),
        "program": lambda pr: (pr[1].get("best_program_P") if pr[1].get("best_program_P") is not None else big),
        "num": lambda pr: pr[1]["num"],
    }.get(sort_by, lambda pr: pr[1]["num"])
    rows.sort(key=keyf)
    out = [r for (_, r) in rows]
    return out[:limit] if limit else out


def get_program(did, condition, program):
    """Full record for one program incl. every trait's program- and regulator-burden + meta."""
    idx = _prog_index(did, condition)
    program = program if str(program).upper().startswith("P") else f"P{program}"
    if program not in idx:
        raise KeyError(f"{program} not in {did}/{condition}")
    return idx[program]


def get_meta(did, condition, program=None):
    """Cross-source random-effects meta (pooled β, P, FDR, I², LOO, source-stratified)."""
    idx = _prog_index(did, condition)
    if program:
        program = program if str(program).upper().startswith("P") else f"P{program}"
        return {"program": program, "condition": condition, "meta": idx[program].get("meta")}
    return [{"program": p["program"], "annotation": p["annotation"], "meta": p.get("meta")}
            for p in _disease(did)["programs"][condition]]


def compare_conditions(did, program):
    """Same program id in rest vs stim48hr. cNMF is INDEPENDENT per condition, so this compares
    two different programs that share a number — the honest caveat is included in the payload."""
    program = program if str(program).upper().startswith("P") else f"P{program}"
    out = {"program": program, "caveat":
           "cNMF is run independently per condition; rest and stim48hr programs with the same number "
           "are NOT the same program. Cross-condition matching is by top-gene Jaccard, not by number.",
           "by_condition": {}}
    for c in _disease(did)["conditions"]:
        idx = _prog_index(did, c)
        if program in idx:
            p = idx[program]
            out["by_condition"][c] = {"annotation": p["annotation"], "top_genes": p["top_genes"][:8],
                                      "meta": p.get("meta")}
    return out


# --------------------------------------------------------------------------- figures / narrative
def list_figures(did, kind=None, trait=None, condition=None):
    figs = _disease(did)["figures"]
    def keep(f):
        return ((kind is None or f["kind"] == kind)
                and (trait is None or f.get("trait") == trait)
                and (condition is None or f.get("condition") == condition))
    return [f for f in figs if keep(f)]


def get_conclusions(did, section="conclusions"):
    return _disease(did)["narrative"].get(section, "")


# --------------------------------------------------------------------------- network / regulators
def list_networks(did):
    nets = _disease(did).get("networks", {})
    return [{"condition": c, "group": g, "label": nets[c][g]["label"],
             "n_programs": nets[c][g]["n_programs"], "n_genes": nets[c][g]["n_genes"]}
            for c in nets for g in nets[c]]


def describe_network(did, condition, group, focus=None):
    """Summary of a (condition, group) tri-partite gene→program→trait network. If `focus` is a program
    or regulator gene, restrict to its neighbourhood. Returns counts + the conserved programs + the
    shared regulators — a compact grounding payload (not the full node list)."""
    net = _disease(did).get("networks", {}).get(condition, {}).get(group)
    if not net:
        raise KeyError(f"no network for {did}/{condition}/{group}")
    nodes = {n["id"]: n for n in net["nodes"]}
    edges = net["edges"]
    if focus:
        keep = {focus}
        for e in edges:
            if e["source"] == focus:
                keep.add(e["target"])
            if e["target"] == focus:
                keep.add(e["source"])
        edges = [e for e in edges if e["source"] in keep and e["target"] in keep]
        nodes = {k: v for k, v in nodes.items() if k in keep}
    progs = sorted([n for n in nodes.values() if n["type"] == "program"],
                   key=lambda n: -n.get("n_traits", 0))
    regs = sorted([n for n in nodes.values() if n["type"] == "gene" and n["role"] == "regulator"],
                  key=lambda n: -n.get("n_traits", 0))
    return {"condition": condition, "group": group, "label": net["label"], "focus": focus,
            "n_programs": sum(1 for n in nodes.values() if n["type"] == "program"),
            "n_genes": sum(1 for n in nodes.values() if n["type"] == "gene"),
            "n_edges": len(edges),
            "conserved_programs": [{"program": p["id"], "annotation": p.get("annotation", ""),
                                    "n_traits": p["n_traits"], "sig": p.get("sig")} for p in progs[:10]],
            "shared_regulators": [{"gene": r["id"], "n_traits": r["n_traits"], "sign": r["sign"]}
                                  for r in regs if r["n_traits"] >= 2][:15]}


def list_regulators(did, shared_only=False, program=None, limit=25):
    regs = _disease(did).get("regulators", [])
    if shared_only:
        regs = [r for r in regs if r.get("shared")]
    if program:
        program = program if str(program).upper().startswith("P") else f"P{program}"
        regs = [r for r in regs if any(t["program"] == program for t in r["targets"])]
    out = [{"gene": r["gene"], "sign": r["sign"], "n_programs": r["n_programs"],
            "n_traits": r["n_traits"], "n_clean_traits": r.get("n_clean_traits"),
            "shared": r.get("shared"), "conditions": r["conditions"],
            "top_targets": [f"{t['program']}({'+' if t['dir'] > 0 else '−'}β{abs(t['beta'] or 0):.1f})"
                            for t in r["targets"][:4]]} for r in regs]
    return out[:limit] if limit else out


def get_regulator(did, gene):
    for r in _disease(did).get("regulators", []):
        if r["gene"].upper() == str(gene).upper():
            return r
    raise KeyError(f"{gene} is not a regulator in {did}")


def trait_panels(did):
    """Headline per (trait, condition): Fig 4C Bonferroni hits + Fig 5a permutation P + QC flag."""
    out = []
    for t in _disease(did)["traits"]:
        out.append({"trait": t["short"], "id": t["id"], "source": t["source"], "flag": t["flag"],
                    "desc": t["desc"], "panels": t["panels"]})
    return out


# --------------------------------------------------------------------------- guardrails (shared)
def system_brief():
    """The non-obvious, load-bearing facts any assistant answering over these results MUST respect.
    Sourced from PLoTM/docs/AGENT_GUIDE.md. Injected into the chat system prompt AND surfaced to
    the MCP client so neither can silently mislead."""
    return (
        "GROUNDING RULES for this dataset (RA_CD4T — patient CD4+ T-cell cNMF, K=30, projected onto "
        "CRISPRi Perturb-seq; do not violate):\n"
        "1. LoF gamma is human population genetics and is CONDITION-INDEPENDENT: one value per gene per "
        "trait, reused for rest and stim48hr. There is no Delta-gamma between conditions. Condition "
        "effects live ONLY in the programs' regulatory betas and burden.\n"
        "2. The programs are ONE patient-derived cNMF decomposition PROJECTED onto Perturb-seq in each "
        "state, so a program has the SAME identity in rest and stim48hr (rest-P6 == stim-P6). Compare "
        "conditions via the per-state regulatory beta / burden of the same program — NOT via re-matching.\n"
        "3. There is a SINGLE genetics source here: GeneBass, 3 traits — GB_RA_custom and GB_RA_M06 "
        "(RA case definitions) and GB_RF (rheumatoid factor, an endophenotype). These are NOT three "
        "independent cohorts (same UK Biobank exomes, overlapping definitions); the 3-trait meta is "
        "random-effects with I^2 shown. There is no Backman/cross-source comparison in this dataset.\n"
        "4. For Fig 5a report the PERMUTATION P, not the Fisher P (programs/regulators are selected on "
        "the same gamma vector, so Fisher is anti-conservative). Fig-5 maps exist for GB_RA_custom and "
        "GB_RA_M06 only; GB_RF has burden/Fig-4C but no Fig-5 network.\n"
        "5. Programs flagged not-estimable (dead / too few control profiles) have no regulator-burden or "
        "meta — exclude them from regulator claims; their program-burden P may still be reported.\n"
        "6. The LoF threshold is a fixed |gamma| = 0.03 (single source, so no cross-source rescaling).\n"
        "7. Programs carry biological annotations (e.g. P6 AP-1/TNFa-NF-kB, P23 Th17/IL-6-STAT3, P2 "
        "Naive-Tcm OXPHOS+translation). Programs of kind TECHNICAL/unresolved are low-confidence — say so."
    )


if __name__ == "__main__":
    # tiny self-check against the first (default) disease
    ds = list_diseases()
    print("diseases:", [d["id"] for d in ds])
    did = ds[0]["id"]
    top = query_programs(did, "rest", sort_by="program", limit=6)
    for r in top:
        fdr = r["meta_FDR"]
        print(f"  {r['program']:>4} progP={r['best_program_P']} metaFDR={fdr if fdr is None else round(fdr, 3)}  {r['annotation']}")
