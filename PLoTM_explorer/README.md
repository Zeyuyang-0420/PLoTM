# PLoTM explorer

An interactive, browser-based way to explore the PLoTM results and **ask questions about them in
natural language**. A Claude chat is grounded in the pipeline's own numbers (via tool calls over the
result tables) *and* in the published literature/databases (via Anthropic server-side web search), so it
can say both "our meta-analysis puts P45 Myeloid activation at FDR 0.028" and "S100A9 is a known RA
synovial marker" — and keep the two clearly separated.

This is a **read-only viewer** for results already produced by [`../PLoTM/`](../PLoTM/); it does not
run the pipeline. It lives in its own directory and touches nothing in `PLoTM/` or `pipeline_results/`.

```
Browser dashboard (program table · figure gallery · chat)
        │  POST /api/chat            GET /api/data
        ▼
FastAPI  server.py   (holds ANTHROPIC_API_KEY — never sent to the browser)
   ├─ local grounding tools ─► proglof_data.py ─► data/demo_data.json
   │     query_programs · get_program · get_meta · compare_conditions
   │     trait_panels · list_figures · get_conclusions
   ├─ web_search (Anthropic server tool) ─► publications + external databases
   └─ Claude Opus tool-use loop ─► grounded answer + tool trace + cited figures
```

## Quick start

```bash
cd PLoTM_explorer

# 1. (already done for RA — regenerate only if pipeline_results/ changed)
PYTHONNOUSERSITE=1 python build_demo_data.py --disease RA

# 2. install + key
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env        # paste your ANTHROPIC_API_KEY

# 3. serve
PYTHONNOUSERSITE=1 uvicorn server:app --port 8080
# open http://localhost:8080
```

The dashboard works without a key (table, figures, conclusions all render from `data/demo_data.json`);
only the **chat** needs `ANTHROPIC_API_KEY`.

## What you can do

- **Programs tab** — all 60 cNMF programs for the chosen condition, sortable by cross-source
  meta-analysis FDR / pooled β / I², filterable to annotated-only or meta-significant, or pivoted to any
  single phenotype's program-burden P. Click a row for a drawer with top genes, the full meta block, and
  the per-phenotype burden table.
- **Network tab** — an interactive version of the combined Fig 5a: a gene→program graph per condition ×
  trait-group (Group1 RA / Group2 RF). Focus on any program or regulator gene to see only its
  neighbourhood; toggle *annotated programs only* or *significant / conserved only*; node size/colour =
  cross-phenotype conservation, ring = γ-sign (genes) or meta-significance (programs), edge width = |β|,
  dashed = program-membership. Click a program → its drawer; click a gene → focus it. Drag nodes to
  declutter.
- **Regulators tab** — every regulator gene (a Stage-B `regulates` edge), ranked by cross-phenotype
  recurrence: γ-sign, #programs, clean/all phenotype counts, and the programs it regulates (β,
  direction). "Shared" surfaces the conserved regulators (e.g. DDX39A, CNIH4, NFIA). Click a row for the
  full target list; jump to the network focused on that gene.
- **Traits tab** — the 8 RA LoF phenotypes × 2 conditions with their Fig 4C hit, **Fig 5a permutation P**
  (not Fisher), and QC flag.
- **Chat** — ask in plain language; answers cite the tools they used (incl. `describe_network`,
  `list_regulators`, `get_regulator`) and any relevant figure, and program references (e.g. `P45`) are
  clickable and cross-highlight the table.

## Grounding & guardrails

The chat's system prompt hard-codes the load-bearing scientific facts (from
[`../PLoTM/docs/AGENT_GUIDE.md`](../PLoTM/docs/AGENT_GUIDE.md), see `proglof_data.system_brief()`):
γ is condition-independent; GeneBass and Backman are the same UK Biobank cohort (not independent); never
compare |γ| across sources with a fixed cutoff; report the permutation P; and honor the QC flags
(`BORDERLINE` / `PRIOR_DOMINATED` / `SIGN_DEGENERATE`). It is told to read the grounding tools before any
quantitative claim.

## Files

| file | role |
|---|---|
| `build_demo_data.py` | flattens `pipeline_results/` → `data/demo_data.json` + copies figures. Disease-keyed. |
| `proglof_data.py` | dependency-free read accessor + shared grounding brief (reused by the Claude Science MCP). |
| `server.py` | FastAPI: static dashboard, `/api/data`, `/api/chat` (Claude tool loop + web search). |
| `static/` | `index.html`, `app.js`, `style.css` — the single-page dashboard (no CDN). |
| `data/demo_data.json` | generated; the entire grounding corpus (~0.4 MB for RA). |

## Adding another disease (SLE-ready)

The schema and UI are disease-keyed. To add a disease, add an entry to `DISEASES` in
`build_demo_data.py` (results root + cNMF spectra paths + phenotype descriptions), run
`python build_demo_data.py --disease <ID>` (or `--all`), and it appears in the disease switcher. No
server or UI change.
