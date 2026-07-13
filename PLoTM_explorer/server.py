"""server.py — PLoTM explorer backend.

Serves the static dashboard and a `/api/chat` endpoint that runs an agentic Claude tool-use loop
grounded in the local pipeline results (via proglof_data.py) plus Anthropic server-side web search for
published knowledge. The ANTHROPIC_API_KEY lives here, server-side — never in the browser.

    pip install -r requirements.txt
    cp .env.example .env   # put your key in it
    PYTHONNOUSERSITE=1 uvicorn server:app --port 8080
    # open http://localhost:8080

Data must be built first:  python build_demo_data.py --disease RA
"""
import json
import os

from fastapi import FastAPI
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

try:
    from dotenv import load_dotenv
    load_dotenv(os.path.join(os.path.dirname(__file__), ".env"))
except Exception:
    pass

import proglof_data as pdata

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL = os.environ.get("PROGLOF_MODEL", "claude-sonnet-5")
MAX_TURNS = 8  # cap the tool loop

app = FastAPI(title="PLoTM explorer")

# ---------------------------------------------------------------------------
# tool definitions exposed to the model (disease is bound per-request, not a tool arg)
# ---------------------------------------------------------------------------
LOCAL_TOOLS = [
    {"name": "query_programs",
     "description": "List/rank cNMF gene programs for a condition. Use to answer 'which programs are "
                    "significant/robust/annotated'. Returns compact rows with meta-analysis stats.",
     "input_schema": {"type": "object", "properties": {
         "condition": {"type": "string", "enum": ["rest", "stim48hr"]},
         "trait": {"type": "string", "description": "optional trait short id, e.g. GB_RF, BK_curated"},
         "annotated_only": {"type": "boolean"},
         "significant": {"type": "string", "enum": ["meta", "program"],
                         "description": "'meta'=cross-source FDR<0.05; 'program'=any trait program-burden P<0.05"},
         "sort_by": {"type": "string", "enum": ["meta_FDR", "meta_pooled", "program", "num"]},
         "limit": {"type": "integer"}},
         "required": ["condition"]}},
    {"name": "get_program",
     "description": "Full record for one program in one condition: annotation, top-loaded genes, "
                    "per-trait program- and regulator-burden, and cross-source meta.",
     "input_schema": {"type": "object", "properties": {
         "condition": {"type": "string", "enum": ["rest", "stim48hr"]},
         "program": {"type": "string", "description": "e.g. P2"}},
         "required": ["condition", "program"]}},
    {"name": "get_meta",
     "description": "Cross-source random-effects meta-analysis (pooled beta, P, FDR, I^2, leave-one-out, "
                    "source-stratified). Omit program to get all programs for the condition.",
     "input_schema": {"type": "object", "properties": {
         "condition": {"type": "string", "enum": ["rest", "stim48hr"]},
         "program": {"type": "string"}},
         "required": ["condition"]}},
    {"name": "compare_conditions",
     "description": "Compare a program id between rest and stim48hr (returns the honest caveat that "
                    "cNMF is independent per condition).",
     "input_schema": {"type": "object", "properties": {"program": {"type": "string"}},
                      "required": ["program"]}},
    {"name": "trait_panels",
     "description": "Headline per (trait, condition): Fig 4C Bonferroni hits, Fig 5a permutation P, and "
                    "the QC flag for each of the 8 RA LoF phenotypes.",
     "input_schema": {"type": "object", "properties": {}}},
    {"name": "list_figures",
     "description": "List available result figures (Fig4C, Fig5a, condition_comparison, combined). "
                    "Return figure file paths the UI can display; cite them in your answer.",
     "input_schema": {"type": "object", "properties": {
         "kind": {"type": "string", "enum": ["Fig4C", "Fig5a", "condition_comparison", "combined"]},
         "trait": {"type": "string"},
         "condition": {"type": "string", "enum": ["rest", "stim48hr"]}}}},
    {"name": "get_conclusions",
     "description": "Return the written conclusions / data-scope narrative for the disease.",
     "input_schema": {"type": "object", "properties": {
         "section": {"type": "string", "enum": ["conclusions", "scope", "overview"]}}}},
    {"name": "describe_network",
     "description": "Summarise the interactive gene→program→trait network for a condition and network "
                    "group. For RA_CD4T the groups are per-trait regulator models (GB_RA_custom, "
                    "GB_RA_M06): the regsubsets regulator programs (with their joint-model weight w_P to "
                    "the trait) + their top trans-regulator genes. Pass `focus` (a program like P18 or a "
                    "regulator gene) to restrict to its neighbourhood.",
     "input_schema": {"type": "object", "properties": {
         "condition": {"type": "string", "enum": ["rest", "stim48hr"]},
         "group": {"type": "string", "description": "network group / trait model id, e.g. GB_RA_custom"},
         "focus": {"type": "string"}}, "required": ["condition", "group"]}},
    {"name": "list_regulators",
     "description": "List regulator genes (genes whose Stage-B knockdown moves a program). Each shows the "
                    "programs it regulates (β, direction) and cross-phenotype recurrence. shared_only=true "
                    "keeps genes recurring across QC-clean phenotypes; filter by program.",
     "input_schema": {"type": "object", "properties": {
         "shared_only": {"type": "boolean"}, "program": {"type": "string"}, "limit": {"type": "integer"}}}},
    {"name": "get_regulator",
     "description": "Full record for one regulator gene: predicted γ-sign and every program it regulates "
                    "with β, direction, trait, condition, and QC flag.",
     "input_schema": {"type": "object", "properties": {"gene": {"type": "string"}},
                      "required": ["gene"]}},
]

WEB_SEARCH_TOOL = {"type": "web_search_20250305", "name": "web_search", "max_uses": 5}


def dispatch(did, name, args):
    a = args or {}
    if name == "query_programs":
        return pdata.query_programs(did, a["condition"], trait=a.get("trait"),
                                    annotated_only=a.get("annotated_only", False),
                                    significant=a.get("significant"),
                                    sort_by=a.get("sort_by", "meta_FDR"), limit=a.get("limit", 25))
    if name == "get_program":
        return pdata.get_program(did, a["condition"], a["program"])
    if name == "get_meta":
        return pdata.get_meta(did, a["condition"], a.get("program"))
    if name == "compare_conditions":
        return pdata.compare_conditions(did, a["program"])
    if name == "trait_panels":
        return pdata.trait_panels(did)
    if name == "list_figures":
        return pdata.list_figures(did, kind=a.get("kind"), trait=a.get("trait"),
                                  condition=a.get("condition"))
    if name == "get_conclusions":
        return pdata.get_conclusions(did, section=a.get("section", "conclusions"))
    if name == "describe_network":
        return pdata.describe_network(did, a["condition"], a["group"], focus=a.get("focus"))
    if name == "list_regulators":
        return pdata.list_regulators(did, shared_only=a.get("shared_only", False),
                                     program=a.get("program"), limit=a.get("limit", 25))
    if name == "get_regulator":
        return pdata.get_regulator(did, a["gene"])
    return {"error": f"unknown tool {name}"}


SYSTEM = (
    "You are the PLoTM results assistant. You help a researcher explore which CD4+ T-cell "
    "transcriptional programs carry a disease's rare loss-of-function (LoF) genetic burden, and which "
    "regulators drive them (Ota et al. 2026 regulators->programs->traits, applied to rheumatoid "
    "arthritis). Answer ONLY from the grounding tools for anything numeric about these results; use "
    "web_search to bring in external published knowledge (papers, databases) about the genes/programs, "
    "and clearly separate 'our result' from 'the literature'. Prefer get_program / get_meta / "
    "query_programs before making any quantitative claim; use describe_network / list_regulators / "
    "get_regulator for the gene->program->trait network or which genes regulate a program. When a figure "
    "is relevant, call list_figures "
    "and name the returned file so it can be shown. Be concise and precise; give the program number AND "
    "its annotation. Never overstate: respect the QC flags and the grounding rules below.\n\n"
    + pdata.system_brief()
)


class ChatIn(BaseModel):
    messages: list          # [{role, content(str)}]
    disease: str = "RA"


def _to_blocks(content):
    """normalise incoming user/assistant content to the API block shape."""
    if isinstance(content, str):
        return [{"type": "text", "text": content}]
    return content


@app.post("/api/chat")
def chat(inp: ChatIn):
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        return JSONResponse(status_code=500,
                            content={"error": "ANTHROPIC_API_KEY not set (see .env.example)."})
    try:
        from anthropic import Anthropic
    except ImportError:
        return JSONResponse(status_code=500,
                            content={"error": "anthropic SDK not installed (pip install -r requirements.txt)."})

    did = inp.disease
    client = Anthropic(api_key=key)
    messages = [{"role": m["role"], "content": _to_blocks(m["content"])} for m in inp.messages]
    tools = LOCAL_TOOLS + [WEB_SEARCH_TOOL]

    trace, cited_figures = [], []
    answer_parts = []
    try:
        for _ in range(MAX_TURNS):
            resp = client.messages.create(model=MODEL, max_tokens=2048, system=SYSTEM,
                                          tools=tools, messages=messages)
            assistant_content = [b.model_dump() for b in resp.content]
            messages.append({"role": "assistant", "content": assistant_content})

            tool_results = []
            for b in assistant_content:
                bt = b.get("type")
                if bt == "text":
                    answer_parts.append(b["text"])
                elif bt == "server_tool_use":            # web_search invocation
                    trace.append({"tool": "web_search", "input": b.get("input", {})})
                elif bt == "tool_use":                   # our local grounding tools
                    trace.append({"tool": b["name"], "input": b.get("input", {})})
                    try:
                        result = dispatch(did, b["name"], b.get("input", {}))
                    except Exception as e:                # noqa: BLE001 — surface to the model
                        result = {"error": str(e)}
                    if b["name"] == "list_figures" and isinstance(result, list):
                        cited_figures.extend(result)
                    tool_results.append({"type": "tool_result", "tool_use_id": b["id"],
                                         "content": json.dumps(result)})

            if resp.stop_reason == "tool_use":
                messages.append({"role": "user", "content": tool_results})
                continue
            if resp.stop_reason == "pause_turn":         # long server-tool turn: continue
                continue
            break
    except Exception as e:  # noqa: BLE001
        return JSONResponse(status_code=502, content={"error": f"Anthropic API error: {e}"})

    # de-dup cited figures by file
    seen, figs = set(), []
    for f in cited_figures:
        if f.get("file") and f["file"] not in seen:
            seen.add(f["file"]); figs.append(f)
    return {"answer": "\n\n".join(answer_parts).strip(), "trace": trace, "figures": figs}


# ---------------------------------------------------------------------------
# static + data
# ---------------------------------------------------------------------------
@app.get("/api/data")
def data():
    p = os.path.join(HERE, "data", "demo_data.json")
    if not os.path.exists(p):
        return JSONResponse(status_code=404,
                            content={"error": "run build_demo_data.py first"})
    return FileResponse(p, media_type="application/json")


@app.get("/")
def index():
    return FileResponse(os.path.join(HERE, "static", "index.html"))


app.mount("/static", StaticFiles(directory=os.path.join(HERE, "static")), name="static")

# The Network tab embeds the React app (../PLoTM_network) built to dist/. Mount it if built.
NET_DIST = os.path.abspath(os.path.join(HERE, "..", "PLoTM_network", "dist"))
if os.path.isdir(NET_DIST):
    app.mount("/network", StaticFiles(directory=NET_DIST, html=True), name="network")
