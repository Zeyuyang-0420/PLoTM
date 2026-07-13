# proglof-mcp — the PLoTM engine as a Claude Science connector

A small **local-command MCP server** (stdio, no third-party dependencies) that gives Claude Science the
seven operations the PLoTM pipeline needs and that **no public database provides**. Everything else an
autonomous run needs — disease background, genetic evidence, phenotype identifiers, the single-cell
substrate — comes from Claude Science's *Featured connectors*, not from here.

## Tools

| tool | wraps | what it does |
|---|---|---|
| `write_config` | — | writes `config/<id>.yaml` from a structured spec (perturb-seq substrate + LoF phenotype panel + params). The single source of truth. |
| `fetch_genebass_burden` | `src/00_download/fetch_genebass_lof.py` | pull GeneBass exome burden + prepare GeneBayes inputs; can append one discovered phenocode first |
| `fetch_backman_burden` | `src/00_download/fetch_backman_lof.py` | pull Backman 2021 burden (via GWAS Catalog) + prepare inputs |
| `run_genebayes` | `run_disease.sh … posterior` | fit gene-level LoF γ + QC gates (GPU) |
| `run_phase` | `drivers/run_disease.sh` | run any pipeline phase: `download·posterior·input·cnmf·stageB·thresholds·burden·figures·annotate·analysis·all` |
| `read_results` | results tables | read `tidy·meta·thresholds·summary` (trait/condition filter + limit) |
| `qc_report` | `lof_thresholds.tsv` | per-trait QC flags + interpretation + the cross-source caveats |

`write_config`'s output round-trips through the pipeline's own `src/cfg_get.py` (verified).

## Register in Claude Science

**Settings → Connectors → Add connector → Local command.**

- **Name:** `proglof`
- **Command:** `python`
- **Arguments (Advanced):** the absolute path to `connectors/proglof-mcp/server.py`
- **Environment (Advanced):** set the variables below.

Then open the connector's **Tools** and set `write_config`, `read_results`, `qc_report`, `run_phase` to
**Always allow** (they are safe and frequently called). Leave `fetch_*`, `run_genebayes` on **Ask each
time** unless you enable *Skip approvals* for the whole connector (only do so once you trust the setup).

`mcp.json` in this folder is the same configuration for MCP clients that read a JSON config (e.g. Claude
Code) — fill in the absolute paths.

## Environment

| var | meaning | default |
|---|---|---|
| `PROGLOF_HOME` | the PLoTM repo | `../../../PLoTM` relative to `server.py` |
| `PROGLOF_PY` | python for the cNMF/scanpy stages (used as `PY` in `run_disease.sh`) | `python` |
| `PROGLOF_RSCRIPT` | Rscript for the R stages | `Rscript` |
| `PROGLOF_RESULTS` | results root for `read_results`/`qc_report` | `$PROGLOF_HOME/results/<id>`, RA falls back to `../pipeline_results` |
| `PROGLOF_REMOTE` | ssh host for heavy phases (cNMF, GeneBayes, Stage B/C) | unset → run locally |
| `PROGLOF_DRY_RUN` | `1` → return the command that *would* run, don't execute | unset |

**Remote compute.** cNMF and GeneBayes need the HPC/GPU box. Two ways to route there: (a) set
`PROGLOF_REMOTE=user@host` so this connector wraps heavy phases in `ssh`; or (b) leave it unset and let
Claude Science's own **remote compute cluster** feature run the jobs (preferred if you configured an SSH
cluster in Settings). Use one, not both.

## Test without Claude Science

```bash
export PROGLOF_HOME=/abs/PLoTM PROGLOF_DRY_RUN=1
printf '%s\n' \
  '{"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}' \
  '{"jsonrpc":"2.0","id":2,"method":"tools/list"}' \
  '{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"qc_report","arguments":{"disease":"RA"}}}' \
  | python server.py
```

You should see the `initialize` handshake, the 7 tools, and the RA QC table.
