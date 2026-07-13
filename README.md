# PLoTM — an interpretable Program · Loss-of-function · Traits Map

**Which transcriptional programs of a cell type carry a disease's rare-variant loss-of-function (LoF)
burden, and which upstream regulators drive them?**

PLoTM operationalises the Ota et al. (2026) *regulators → programs → traits* model
([10.1038/s41586-025-09866-3](https://doi.org/10.1038/s41586-025-09866-3)) as a reproducible,
config-driven, agent-runnable pipeline. It joins a CRISPRi **Perturb-seq** atlas (cNMF gene programs +
per-gene regulatory effects) to **gene-level LoF posteriors** from two independent exome-burden sources
(GeneBass and Backman 2021), and asks — with permutation nulls, FDR, and cross-source meta-analysis —
which programs are enriched for a disease's genetic risk. The reference configuration reproduces a
**rheumatoid-arthritis (RA)** analysis on primary human **CD4⁺ T cells**.

## The four parts

| directory | what it is |
|---|---|
| [PLoTM/](PLoTM/README.md) | **The compute engine.** The config-driven `regulators → programs → traits` pipeline (cNMF · Stage B/C burden · figures · meta-analysis). Start here. |
| [PLoTM_agent/](PLoTM_agent/README.md) | **Autonomous disease research on Claude Science.** A skill + one MCP connector that turn Claude into an agent which researches a *new* disease end-to-end from a single prompt. |
| [PLoTM_explorer/](PLoTM_explorer/README.md) | **Natural-language results explorer.** A FastAPI + Claude dashboard (program table · figure gallery · chat) grounded in the pipeline's own numbers and the literature. |
| [PLoTM_network/](PLoTM_network/README.md) | **Interactive gene → program → trait graph.** A React + D3-force network view of the results; self-contained, no server. |

## Quick start

The engine's reference RA run and per-stage commands are documented in **[PLoTM/README.md](PLoTM/README.md)**.
Each app directory has its own README with setup and run instructions.

## Credits & citation

Method: Ota et al., *Causal modelling of gene effects from regulators to programs to traits*, Nature
(2026), doi:10.1038/s41586-025-09866-3. Perturb-seq: CZI genome-scale CD4⁺ T-cell Perturb-seq
(Zhu, Dann, Yan et al.; Pritchard/Marson labs). LoF: GeneBass (Karczewski et al. 2022) and Backman
et al. 2021 (via GWAS Catalog). cNMF: Kotliar et al. 2019. GeneBayes: Zeng et al.
