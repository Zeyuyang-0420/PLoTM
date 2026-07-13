# PLoTM on Claude Science — autonomous disease LoF research

Turn **Claude Science** (Anthropic's desktop science workbench) into an agent that researches a new
disease end-to-end with the PLoTM pipeline: it reads the disease's background from the literature,
discovers the relevant LoF phenotypes, runs the whole *regulators → programs → traits* analysis, and
writes an interpreted report — from a single prompt.

The insight is that **Claude Science already provides the agent, the UI, the sandbox, the literature
access, and the scientific-database connectors.** So this package is small: a **skill** that encodes the
decision procedure and guardrails, **one custom connector** that exposes the PLoTM engine, and the
setup that wires them to compute.

```
Claude Science (desktop app = agent + UI + sandbox + reviewer)
   ├─ Featured connectors ──► discovery, no code from us:
   │     OLS/MyGene (disease→ontology→genes) · GWAS Catalog/FinnGen/Open Targets (genetics + phenotype IDs)
   │     literature access + OpenAlex/PubMed/bioRxiv (background) · CELLxGENE/GEO (single-cell substrate)
   ├─ skill: disease-lof-research ──► the 6-step loop + the load-bearing guardrails
   ├─ connector: proglof-mcp ──► write_config · fetch_{genebass,backman}_burden · run_genebayes
   │                              · run_phase · read_results · qc_report   (wraps ../PLoTM)
   └─ SSH remote compute ──► cNMF · GeneBayes · Stage B/C on the HPC/GPU box
```

## Contents

| path | what |
|---|---|
| `skills/disease-lof-research/` | the Agent Skill: `SKILL.md` + `reference/{connectors,guardrails,pipeline,report_template}.md` |
| `connectors/proglof-mcp/` | the PLoTM-engine MCP server (`server.py`, no deps) + `README.md` + `mcp.json` |
| `claude-science-setup.md` | one-time setup: grant folders · enable Featured connectors · register proglof-mcp · attach SSH compute · add the skill |
| `project/research_brief_template.md` | the working record the skill maintains per disease |
| `prompts/kickoff.md` | the single prompt that starts an autonomous run |
| `examples/T1D/` | a worked reference (type 1 diabetes) + the config it should produce — the acceptance target |

## Quick start

1. Follow `claude-science-setup.md` (once).
2. Open a session in the project and paste `prompts/kickoff.md` with your disease.
3. Claude runs the skill: background → phenotype panel → `write_config` → pipeline via `proglof-mcp` →
   `report.md`. It asks you to approve each new folder/host/job; the reviewer then checks the claims.

## Why the guardrails matter

The skill hard-codes the facts that have burned real runs (see `skills/…/reference/guardrails.md`): γ is
condition-independent; GeneBass and Backman are one UK Biobank cohort (not independent replication); never
share a fixed |γ| cutoff across sources; report the permutation P not Fisher; honor the QC gates. Without
them an agent will confidently over-claim. With them, an underpowered rare-LoF disease (like SLE's IFN
axis, or likely T1D) is reported honestly as null rather than dressed up.

## Relationship to the rest of the repo

- `../PLoTM/` — the pipeline this drives (config-driven `run_disease.sh`, the stages, `cfg_get.py`).
- `../pipeline_results/` — the RA reference results (`read_results` reads these for RA).
- `../PLoTM_explorer/` — the *other* deliverable: an interactive web explorer for results already
  produced. This package instead *produces* results for a new disease.

Requires the Claude Science desktop app (beta, macOS/Linux, Pro/Max/Team/Enterprise).
