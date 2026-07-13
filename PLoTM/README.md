# PLoTM: An Interpretable Program–LoF–Trait Map for Discovering Disease Regulatory Programs and Regulators in Rheumatoid Arthritis

**Which transcriptional programs of a cell type carry a disease's rare-variant loss-of-function (LoF)
burden, and which upstream regulators drive them?**

PLoTM (Programs Loss-of-function Traits Map) operationalises the Ota et al. (2026) *regulators → programs → traits* model
([10.1038/s41586-025-09866-3](https://doi.org/10.1038/s41586-025-09866-3)) as a reproducible,
**config-driven, agent-runnable** pipeline. It joins a CRISPRi **Perturb-seq** atlas (cNMF gene
programs + per-gene regulatory effects) to **gene-level LoF posteriors** for a panel of disease
phenotypes from **two independent exome-burden sources** (GeneBass and Backman 2021), and asks — with
proper permutation nulls, FDR, and cross-source meta-analysis — which programs are enriched for the
disease's genetic risk.

The reference configuration reproduces a **rheumatoid-arthritis (RA)** analysis on primary human
**CD4⁺ T cells** across two stimulation conditions (Rest, Stim48hr).

**Two modes** (the second is optional):

| mode | programs learned from | tag | config |
|---|---|---|---|
| **pure Perturb-seq** (default) | the Perturb-seq pseudobulk itself | `proglof-v1.0-pure-perturbseq` | `config/RA.yaml` |
| **disease single-cell basis** (optional) | disease single-cell data (e.g. patient CD4⁺ T cells); Perturb-seq is *projected* onto that fixed basis | `proglof-v1.1-disease-singlecell` | `config/RA_singlecell.yaml` |

The second mode is switched on purely by adding a `disease_singlecell:` block to the config and
shifts the question from *"Perturb programs → trait"* to *"Perturb regulators → **disease** programs
→ trait"*. See [docs/DISEASE_BASIS.md](docs/DISEASE_BASIS.md).

---

## The idea in one diagram

**Where the gene programs come from** — pick one basis, everything downstream is identical:

```
  A ─ pure Perturb-seq basis (default, v1.0)
      pseudobulk h5ad ─► cNMF (K=60) ─► gene programs ──┐
                                        └► usages ──────┤

  B ─ disease single-cell basis (optional, v1.1)
      disease scRNA (e.g. RA CD4⁺ T) ─► cNMF (K=30) ─► FIXED gene programs (basis W)
              │                                              │
              └──────── project Perturb-seq pseudobulk onto W (NNLS, U≥0) ─► new usages ──┐
                        (program identities = the disease's; usages = per condition)      │
```

Both feed the **same** disease-independent core, then join the LoF genetics:

```
                 programs + usages (from A or B, per condition)       LoF genetics (disease, 2 sources)
                 ─────────────────────────────────────────────       ─────────────────────────────────
   gene programs ─────────────────────────────────┐        exome burden ─► GeneBayes ─► gamma (per gene)
   usages ─► Stage B ─► beta_x(P)  (regulatory effects)     │                             │
                                                   │        │                             │
                                (LoF-independent, computed ONCE per condition)            │
                                                   ▼                                       ▼
      Stage C burden  ──►  Figure 1   (program-burden vs regulator-burden)  ◄── per-trait gamma
      Fig-2 model     ──►  Figure 2   (regulators → programs → trait + permutation P)
      annotation      ──►  program labels (GO / hallmark / markers)
      condition + meta ─►  Rest-vs-Stim comparison, cross-source meta-analysis, conserved-gene networks
```

The cNMF programs and Stage-B β are **independent of the disease's genetics** — computed once per
condition, then every LoF phenotype reuses the identical matrix (enforced by symlinks). Swapping the
disease is swapping the phenotype list. In mode **B**, the programs additionally carry the disease's
own single-cell structure: the question shifts from *"Perturb programs → trait"* to *"Perturb
regulators → **disease** programs → trait"*. See [docs/DISEASE_BASIS.md](docs/DISEASE_BASIS.md).

---

## Repository layout

```
config/            disease.example.yaml  — the single file an agent edits per disease
drivers/           run_disease.sh (orchestrator) · run_batch.sh · run_trait.sh · compute_thresholds.py
src/
  00_download/     perturb-seq + LoF burden from GeneBass and Backman (2 sources)
  01_lof_posterior/GeneBayes gamma (GPU) + QC gates + multi-run ensembling
  02_input_prep/   pseudobulk -> cNMF input h5ad + covariate metadata
  03_cnmf_stageA/  cNMF program discovery (paper-faithful)
  04_regulatory_stageB/  beta_x(P) regulatory effect sizes
  05_burden_stageC/      program & regulator burden tests (Fig-1 inputs)
  06_annotation/   program annotation + curated labels
  07_figures/      Fig 1, Fig 2 (2 steps), multi-trait networks, summaries
  08_condition_analysis/ tidy table, meta-analysis, condition interaction, permutation, figures
environment/       cnmf.yml · r_core.yml · genebayes.yml
docs/              PIPELINE.md · AGENT_GUIDE.md
```

---

## Quick start (reference RA run)

```bash
# 0. environments
conda env create -p ./envs/cnmf      -f environment/cnmf.yml
conda env create -p ./envs/r_core    -f environment/r_core.yml
conda env create -p ./envs/genebayes -f environment/genebayes.yml   # needs a CUDA GPU

# 1. point the config at your data, then run the whole thing
export PY=./envs/cnmf/bin/python  RSCRIPT=./envs/r_core/bin/Rscript
drivers/run_disease.sh config/RA.yaml all

# …or run one stage at a time (download | posterior | input | cnmf | stageB |
#    thresholds | burden | figures | annotate | analysis)
drivers/run_disease.sh config/RA.yaml cnmf
```

Every knob (dataset URLs, phenotype list, K, thresholds, compute) lives in the YAML — see
`config/disease.example.yaml`.

---

## Running a NEW disease (agent-based)

The pipeline is built to be driven by an agent primarily from a new config. In short: copy
`config/disease.example.yaml`, set the perturb-seq dataset + `obs_map`, list the disease phenotypes and
their GeneBass/Backman identifiers, pick the conditions, and call `run_disease.sh <config> all`. The
upstream stages are fully config-driven; the reference burden/figure/analysis stages run as-is for RA
and take new paths via env vars; only the two `download`/`posterior` fetch scripts need a per-disease
edit to supply phenotype identifiers (see the INTEGRATION STATUS block in `drivers/run_disease.sh`). QC
gates flag unreliable posteriors automatically; the driver is resumable and skips completed work. The
full decision procedure — including the non-obvious data facts an agent must respect — is in
**[docs/AGENT_GUIDE.md](docs/AGENT_GUIDE.md)**.

---

## What the reference RA run found (headline)

* **Robust cross-source programs** (random-effects meta across 5 RA phenotypes, I²≈0, FDR-significant,
  leave-one-phenotype-out stable, GeneBass and Backman agree): **Translation initiation, TNF-NFkB,
  Myeloid activation (S100A9⁺), IFN-α** — canonical RA biology.
* **Shared regulator genes exceed chance** (permutation P < 1e-3): **ATXN1, NFIA** recur across
  phenotypes and conditions. Shared *programs*, by contrast, do **not** exceed chance.
* **Condition dependence is weak**: LoF γ is population genetics and is condition-independent, so
  associations are mostly constitutive; the one Rest-specific effect is Prostaglandin-RANKL-Th2.
* Two Backman sub-phenotypes are auto-flagged (prior-dominated / sign-degenerate) and excluded from
  meta-analysis; a putative "olfactory/noise" program is QC-flagged and excluded from interpretation.

Full write-up: `docs/` and the reference results under `results/RA/` (produced by a run).

---

## The rest of PLoTM (sibling projects)

This repo (`PLoTM/`) is the **compute engine**. Three companion projects wrap it into an end-to-end
system — from an autonomous agent that runs a new disease, to interactive ways of reading the results:

| project | what it is |
|---|---|
| [../PLoTM_agent/](../PLoTM_agent/README.md) | **Autonomous disease research on Claude Science.** A skill + one MCP connector that turn Claude into an agent which researches a *new* disease end-to-end from a single prompt: literature background → LoF phenotype discovery → `write_config` → runs this pipeline → interpreted report. |
| [../PLoTM_explorer/](../PLoTM_explorer/README.md) | **Natural-language results explorer.** A FastAPI + Claude dashboard (program table · figure gallery · chat) grounded in the pipeline's own numbers *and* the literature via web search — a read-only viewer over `pipeline_results/`. |
| [../PLoTM_network/](../PLoTM_network/README.md) | **Interactive gene → program → trait graph.** A React + D3-force network view that replaces the static Fig-2 maps with a physics-based, draggable, zoomable three-layer graph; fully self-contained, no server. |

---

## Documentation

| doc | what |
|---|---|
| [docs/PIPELINE.md](docs/PIPELINE.md) | stage-by-stage methodology, inputs/outputs, parameters |
| [docs/AGENT_GUIDE.md](docs/AGENT_GUIDE.md) | how to adapt to a new disease; data facts an agent must respect |

## Credits & citation

Method: Ota et al., *Causal modelling of gene effects from regulators to programs to traits*, Nature
(2026), doi:10.1038/s41586-025-09866-3. Perturb-seq: CZI genome-scale CD4⁺ T-cell Perturb-seq
(Zhu, Dann, Yan et al.; Pritchard/Marson labs). LoF: GeneBass (Karczewski et al. 2022) and Backman
et al. 2021 (via GWAS Catalog). cNMF: Kotliar et al. 2019. GeneBayes: Zeng et al.
