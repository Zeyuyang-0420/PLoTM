#!/usr/bin/env python3
"""proglof-mcp — a Model Context Protocol server that exposes the PLoTM engine to Claude Science.

Claude Science's *Featured connectors* already cover discovery (GWAS Catalog / Open Targets / OLS /
literature / CELLxGENE …). This connector adds only the operations that no public database provides:
writing a disease config, fetching the LoF burden, fitting GeneBayes, running the pipeline phases, and
reading the results back. Register it in Claude Science: Settings > Connectors > Add > Local command,
command = `python /abs/path/to/server.py`.

Transport: newline-delimited JSON-RPC 2.0 over stdio (the MCP stdio transport). No third-party deps, so
it runs in the CS sandbox as-is and is unit-testable by piping.

Environment (set these in the connector's Advanced > Environment):
  PROGLOF_HOME     path to the PLoTM repo               (default: ../../../PLoTM from here)
  PROGLOF_PY       python for run_disease.sh (scanpy/cnmf)(default: `python`)
  PROGLOF_RSCRIPT  Rscript for the R stages               (default: `Rscript`)
  PROGLOF_RESULTS  results root for the read tools         (default: $PROGLOF_HOME/results, with an
                                                            RA fallback to <repo>/../pipeline_results)
  PROGLOF_REMOTE   optional ssh host; if set, heavy phases run there (Claude Science SSH cluster)
  PROGLOF_DRY_RUN  "1" => return the command that WOULD run instead of executing (safe preview/tests)
"""
import json
import os
import subprocess
import sys

# ---------------------------------------------------------------------------
_SELF = os.path.dirname(os.path.abspath(__file__))
HOME = os.environ.get("PROGLOF_HOME", os.path.abspath(os.path.join(_SELF, "..", "..", "..", "PLoTM")))
PY = os.environ.get("PROGLOF_PY", "python")
RSCRIPT = os.environ.get("PROGLOF_RSCRIPT", "Rscript")
REMOTE = os.environ.get("PROGLOF_REMOTE", "")
DRY = os.environ.get("PROGLOF_DRY_RUN", "") in ("1", "true", "yes")
PROTOCOL = "2024-11-05"
VERSION = "0.1.0"


def results_root(disease):
    if os.environ.get("PROGLOF_RESULTS"):
        return os.environ["PROGLOF_RESULTS"]
    cand = os.path.join(HOME, "results", disease)
    if os.path.isdir(cand):
        return cand
    # RA reference results live beside the repo
    ra = os.path.abspath(os.path.join(HOME, "..", "pipeline_results"))
    return ra if (disease.upper() == "RA" and os.path.isdir(ra)) else cand


# --------------------------------------------------------------------------- YAML emit
def _yscalar(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    s = str(v)
    return '"%s"' % s if (s == "" or any(c in s for c in ' :#{}[],&*?|<>=!%@`"\'') ) else s


def emit_yaml(obj, indent=0):
    """Minimal YAML emitter for the PLoTM config shape (nested maps, scalar lists inline,
    dict lists as inline-dict block items). Round-trips through cfg_get.py."""
    pad = "  " * indent
    out = []
    for k, v in obj.items():
        if isinstance(v, dict):
            out.append(f"{pad}{k}:")
            out.append(emit_yaml(v, indent + 1))
        elif isinstance(v, list):
            if v and isinstance(v[0], dict):
                out.append(f"{pad}{k}:")
                for item in v:
                    inner = ", ".join(f"{ik}: {_yscalar(iv)}" for ik, iv in item.items())
                    out.append(f"{pad}  - {{ {inner} }}")
            else:
                out.append(f"{pad}{k}: [{', '.join(_yscalar(x) for x in v)}]")
        else:
            out.append(f"{pad}{k}: {_yscalar(v)}")
    return "\n".join(x for x in out if x != "")


# --------------------------------------------------------------------------- shell
def _run(cmd, cwd=None, phase_heavy=False):
    """Execute a shell command, or (DRY) return it. Heavy phases route to PROGLOF_REMOTE if set."""
    if phase_heavy and REMOTE:
        cmd = f"ssh {REMOTE} 'cd {cwd or HOME} && {cmd}'"
        cwd = None
    if DRY:
        return {"dry_run": True, "command": cmd, "cwd": cwd or HOME}
    p = subprocess.run(cmd, shell=True, cwd=cwd or HOME, capture_output=True, text=True)
    tail = lambda s: "\n".join(s.strip().splitlines()[-40:])
    return {"command": cmd, "returncode": p.returncode,
            "stdout_tail": tail(p.stdout), "stderr_tail": tail(p.stderr)}


def _read_tsv(path, limit=None, where=None):
    if not os.path.exists(path):
        return {"error": f"not found: {path}"}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            vals = line.rstrip("\n").split("\t")
            row = dict(zip(header, vals))
            if where and not all(row.get(k) == v for k, v in where.items()):
                continue
            rows.append(row)
            if limit and len(rows) >= limit:
                break
    return {"path": path, "columns": header, "n": len(rows), "rows": rows}


def _driver(config, phase, heavy=False):
    env = f"PY={PY} RSCRIPT={RSCRIPT}"
    cmd = f"{env} bash drivers/run_disease.sh {config} {phase}"
    return _run(cmd, cwd=HOME, phase_heavy=heavy)


# --------------------------------------------------------------------------- tools
def t_write_config(a):
    disease = a["disease"]
    did = disease["id"] if isinstance(disease, dict) else disease
    cfg = {
        "disease": disease if isinstance(disease, dict) else {"id": did, "name": a.get("name", did)},
        "perturbseq": a["perturbseq"],
        "lof": {"posterior_dir": a.get("posterior_dir", "data/lof/GeneBayes_posterior"),
                "traits": a["traits"]},
        "params": a.get("params", {"K": 60, "lof_threshold": "q99", "stageB_min_profiles": 4,
                                    "fig5": {"PN": 5, "RN": 3, "TOP": 200, "NPERM": 5000,
                                             "NPERM_OBS": 10000}}),
        "compute": a.get("compute", {"ncores": 64, "throttle": 80, "gpu": True}),
    }
    if a.get("trait_groups"):
        cfg["params"]["trait_groups"] = a["trait_groups"]
    text = emit_yaml(cfg) + "\n"
    out = os.path.join(HOME, "config", f"{did}.yaml")
    if DRY:
        return {"dry_run": True, "path": out, "yaml": text}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        fh.write(text)
    return {"path": out, "n_traits": len(cfg["lof"]["traits"]),
            "conditions": cfg["perturbseq"].get("conditions"), "yaml": text}


def t_fetch_genebass(a):
    cfg = a["config"]
    if a.get("phenocode") and a.get("name"):
        _append_trait(cfg, {"name": a["name"], "source": "genebass",
                            "phenocode": a["phenocode"], "class": a.get("class", "binary")})
    return {"note": "fetch GeneBass burden for genebass traits in the config",
            **_run(f"{PY} src/00_download/fetch_genebass_lof.py --config {cfg} && "
                   f"{PY} src/00_download/prepare_lof_inputs_genebass.py --config {cfg}", cwd=HOME)}


def t_fetch_backman(a):
    cfg = a["config"]
    if a.get("accession") and a.get("name"):
        _append_trait(cfg, {"name": a["name"], "source": "backman",
                            "accession": a["accession"], "class": a.get("class", "binary")})
    return {"note": "fetch Backman burden (via GWAS Catalog) for backman traits in the config",
            **_run(f"{PY} src/00_download/fetch_backman_lof.py --config {cfg} && "
                   f"{PY} src/00_download/prepare_lof_inputs_backman.py --config {cfg}", cwd=HOME)}


def _append_trait(cfg_path, trait):
    """Append one discovered phenotype to an existing config's lof.traits (idempotent by name)."""
    p = cfg_path if os.path.isabs(cfg_path) else os.path.join(HOME, cfg_path)
    if not os.path.exists(p):
        return
    line = "    - { " + ", ".join(f"{k}: {_yscalar(v)}" for k, v in trait.items()) + " }\n"
    with open(p) as fh:
        txt = fh.read()
    if f"name: {trait['name']}" in txt:
        return
    # insert after the 'traits:' line
    out, done = [], False
    for ln in txt.splitlines(keepends=True):
        out.append(ln)
        if not done and ln.strip().startswith("traits:"):
            out.append(line); done = True
    with open(p, "w") as fh:
        fh.write("".join(out))


def t_run_genebayes(a):
    return {"note": "GeneBayes posterior + QC gates (GPU; routed to PROGLOF_REMOTE if set)",
            **_driver(a["config"], "posterior", heavy=True)}


def t_run_phase(a):
    phase = a["phase"]
    heavy = phase in ("cnmf", "stageB", "burden", "figures", "posterior", "all")
    return _driver(a["config"], phase, heavy=heavy)


def t_read_results(a):
    disease = a.get("disease", "RA")
    kind = a.get("kind", "meta")
    rr = results_root(disease)
    paths = {
        "tidy": os.path.join(rr, "condition_comparison/tables/tidy_results.tsv"),
        "meta": os.path.join(rr, "condition_comparison/tables/meta_regulator_RA5.tsv"),
        "thresholds": os.path.join(rr, "lof_thresholds.tsv"),
        "summary": os.path.join(rr, "summary.tsv"),
    }
    if kind not in paths:
        return {"error": f"kind must be one of {list(paths)}"}
    where = {}
    if a.get("condition"):
        where["condition"] = a["condition"]
    if a.get("trait"):
        where["trait"] = a["trait"]
    return _read_tsv(paths[kind], limit=a.get("limit", 100), where=where or None)


def t_qc_report(a):
    disease = a.get("disease", "RA")
    rr = results_root(disease)
    res = _read_tsv(os.path.join(rr, "lof_thresholds.tsv"))
    if "rows" in res:
        notes = {
            "PRIOR_DOMINATED": "trait data barely move the posterior (corr>0.95) — exclude from robust claims",
            "SIGN_DEGENERATE": "sign concordance collapses (<2% or >98% negative gamma) — Fig 5a not interpretable",
            "BORDERLINE": "near a QC gate — treat with caution",
        }
        res["interpretation"] = notes
        res["reminder"] = ("gamma is condition-independent; GeneBass & Backman are the same UK Biobank "
                           "cohort; never share a fixed |gamma| cutoff across sources (per-trait q99).")
    return res


TOOLS = [
    {"name": "write_config",
     "description": "Write a PLoTM disease config (config/<id>.yaml) from a structured spec: the "
                    "perturb-seq substrate + obs_map, the LoF phenotype panel (GeneBass phenocodes / "
                    "Backman accessions discovered via the Featured connectors), and params. This is "
                    "the single source of truth the pipeline runs from.",
     "inputSchema": {"type": "object", "properties": {
         "disease": {"type": "object", "description": "{id, name}"},
         "perturbseq": {"type": "object", "description": "{url, pseudobulk_h5ad, obs_map, conditions[]}"},
         "traits": {"type": "array", "items": {"type": "object"},
                    "description": "[{name, source: genebass|backman, phenocode|accession, class}]"},
         "params": {"type": "object"}, "trait_groups": {"type": "object"},
         "posterior_dir": {"type": "string"}},
         "required": ["disease", "perturbseq", "traits"]}},
    {"name": "fetch_genebass_burden",
     "description": "Fetch the GeneBass exome gene-burden summary stats for the genebass traits in a "
                    "config and prepare the GeneBayes inputs. Optionally append one discovered "
                    "phenotype (phenocode+name) to the config first.",
     "inputSchema": {"type": "object", "properties": {
         "config": {"type": "string"}, "phenocode": {"type": "string"}, "name": {"type": "string"},
         "class": {"type": "string"}}, "required": ["config"]}},
    {"name": "fetch_backman_burden",
     "description": "Fetch the Backman 2021 exome burden (via GWAS Catalog) for the backman traits in a "
                    "config and prepare the GeneBayes inputs. Optionally append one discovered "
                    "phenotype (accession+name) first.",
     "inputSchema": {"type": "object", "properties": {
         "config": {"type": "string"}, "accession": {"type": "string"}, "name": {"type": "string"},
         "class": {"type": "string"}}, "required": ["config"]}},
    {"name": "run_genebayes",
     "description": "Fit the GeneBayes gene-level LoF posterior (gamma) per trait and run the QC gates. "
                    "GPU; routed to the SSH compute cluster if PROGLOF_REMOTE is set.",
     "inputSchema": {"type": "object", "properties": {"config": {"type": "string"}},
                     "required": ["config"]}},
    {"name": "run_phase",
     "description": "Run one PLoTM pipeline phase via drivers/run_disease.sh. Phases: download | "
                    "posterior | input | cnmf | stageB | thresholds | burden | figures | annotate | "
                    "analysis | all. Heavy phases route to PROGLOF_REMOTE if set. Resumable.",
     "inputSchema": {"type": "object", "properties": {
         "config": {"type": "string"},
         "phase": {"type": "string", "enum": ["download", "posterior", "input", "cnmf", "stageB",
                    "thresholds", "burden", "figures", "annotate", "analysis", "all"]}},
         "required": ["config", "phase"]}},
    {"name": "read_results",
     "description": "Read a PLoTM results table: kind = tidy (per program x trait x condition) | meta "
                    "(cross-source meta-analysis) | thresholds | summary. Optional trait/condition "
                    "filter + limit.",
     "inputSchema": {"type": "object", "properties": {
         "disease": {"type": "string"},
         "kind": {"type": "string", "enum": ["tidy", "meta", "thresholds", "summary"]},
         "trait": {"type": "string"}, "condition": {"type": "string"}, "limit": {"type": "integer"}},
         "required": ["disease", "kind"]}},
    {"name": "qc_report",
     "description": "Return the per-trait QC table (thresholds + PRIOR_DOMINATED / SIGN_DEGENERATE / "
                    "BORDERLINE flags) plus the interpretation and the cross-source caveats an "
                    "interpreter must respect.",
     "inputSchema": {"type": "object", "properties": {"disease": {"type": "string"}},
                     "required": ["disease"]}},
]

HANDLERS = {
    "write_config": t_write_config, "fetch_genebass_burden": t_fetch_genebass,
    "fetch_backman_burden": t_fetch_backman, "run_genebayes": t_run_genebayes,
    "run_phase": t_run_phase, "read_results": t_read_results, "qc_report": t_qc_report,
}


def call_tool(name, args):
    if name not in HANDLERS:
        return {"content": [{"type": "text", "text": f"unknown tool {name}"}], "isError": True}
    try:
        result = HANDLERS[name](args or {})
        return {"content": [{"type": "text", "text": json.dumps(result, indent=2)}]}
    except Exception as e:  # noqa: BLE001
        return {"content": [{"type": "text", "text": f"error in {name}: {e}"}], "isError": True}


# --------------------------------------------------------------------------- stdio loop
def _send(obj):
    sys.stdout.write(json.dumps(obj) + "\n")
    sys.stdout.flush()


def serve():
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue
        method, mid, params = msg.get("method"), msg.get("id"), msg.get("params") or {}
        if method == "initialize":
            result = {"protocolVersion": params.get("protocolVersion", PROTOCOL),
                      "capabilities": {"tools": {}},
                      "serverInfo": {"name": "proglof-mcp", "version": VERSION}}
        elif method == "tools/list":
            result = {"tools": TOOLS}
        elif method == "tools/call":
            result = call_tool(params.get("name"), params.get("arguments"))
        elif method == "ping":
            result = {}
        elif method and method.startswith("notifications/"):
            continue  # notifications get no response
        else:
            if mid is not None:
                _send({"jsonrpc": "2.0", "id": mid,
                       "error": {"code": -32601, "message": f"method not found: {method}"}})
            continue
        if mid is not None:
            _send({"jsonrpc": "2.0", "id": mid, "result": result})


if __name__ == "__main__":
    serve()
