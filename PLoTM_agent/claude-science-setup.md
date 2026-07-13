# Setting up Claude Science to auto-research a new disease

One-time setup. After this, a single prompt (`prompts/kickoff.md`) runs the whole pipeline for any disease.

Requires the Claude Science desktop app (macOS 13+/Linux x64, Pro/Max/Team/Enterprise plan). Install and
sign in per <https://claude.com/docs/claude-science/get-started>, then `claude-science serve`.

## 1. Grant the project folders

Create a project (Core concepts → projects group sessions + artifacts). Grant:

- **read** `…/PLoTM/` (the pipeline) and `…/pipeline_results/` (the RA reference results),
- **read+write** a workspace, e.g. `…/PLoTM_research/` — the skill writes `research/<disease>/` here.

You approve each new folder the first time Claude uses it.

## 2. Enable the Featured connectors (Settings → Connectors)

These do all the discovery — turn on (most are on by default):

- **Genes & Ontologies** (OLS → MONDO/EFO, MyGene, GO, Reactome) — disease→ontology→genes, program biology
- **Human Genetics** (GWAS Catalog, FinnGen, BioBank Japan) — genetic architecture + Backman accessions
- **Clinical Genomics** (Open Targets, ClinGen, CIViC) — trait↔gene evidence, phenotype mapping
- **Variants** (gnomAD, ClinVar) — LoF-intolerance / known-gene sanity checks
- **Literature Graph** (OpenAlex, arXiv) + **CellGuide** (CELLxGENE) + **GEO/ArrayExpress** — background + substrate

Add the **directory connectors** **PubMed** and **bioRxiv** (Browse Connectors Directory; on Team/Enterprise
an admin adds them). Set a **contact email** (Settings → General) to enable the Unpaywall literature step.

## 3. Register the pipeline connector (proglof-mcp)

Settings → Connectors → **Add connector → Local command**:

- **Name** `proglof`
- **Command** `python`
- **Arguments** `/abs/path/PLoTM_agent/connectors/proglof-mcp/server.py`
- **Environment** `PROGLOF_HOME`, `PROGLOF_PY`, `PROGLOF_RSCRIPT` (see `connectors/proglof-mcp/README.md`)

On the connector's **Tools** page set `write_config`, `read_results`, `qc_report`, `run_phase` to **Always
allow**. (Environment variables for local connectors are stored unencrypted, readable only by your
account — fine for paths, don't put secrets there.)

## 4. Attach the compute (cNMF + GeneBayes are heavy)

Claude Science → **remote compute clusters**: add the HPC/GPU box as an SSH host and let Claude submit the
`cnmf` / `stageB` / `posterior` phases there. *Either* use this *or* set `PROGLOF_REMOTE=user@host` on the
proglof connector so it wraps heavy phases in ssh itself — not both. GeneBayes needs a CUDA GPU.

## 5. Add the skill

Settings → **Skills → Add skill → Import from GitHub** (or **Upload a skill**) → point at
`PLoTM_agent/skills/disease-lof-research/`. It loads automatically when you ask to research a
disease, or insert it with `/` in the composer. (Import from GitHub needs a token in Settings → Credentials
for a private repo.)

## 6. Run

Open a session in the project and paste `prompts/kickoff.md` with the disease name. Claude follows the
skill: background → phenotype panel → config → pipeline → interpreted report. It asks for approval at each
new folder/host/job (that is by design); the built-in **reviewer** then checks the report's claims against
the execution record.

## Autonomy vs. approvals

Claude Science is approval-gated on purpose. For a near-unattended run, approve the proglof connector's
safe tools as *Always allow* (step 3) and grant folders/hosts at *This project* scope. Do **not** launch
`claude-science serve --dangerously-skip-approvals` for real work — it approves everything, for everything.
