---
name: disease-lof-research
description: >
  Run the PLoTM regulators→programs→traits analysis end-to-end for a disease the user names (e.g.
  "research multiple sclerosis", "which CD4 T-cell programs carry the type-1-diabetes LoF burden").
  Autonomously builds disease background from the literature, selects LoF phenotypes from GeneBass AND
  Backman using the human-genetics connectors, writes a PLoTM config, runs the pipeline through the
  proglof-mcp connector (cNMF → Stage B/C → figures), and writes an interpreted report that separates
  robust cross-source findings from phenotype-specific or QC-flagged ones. Use whenever the goal is to
  apply the CD4 T-cell Perturb-seq × rare-variant-burden pipeline to a new disease or phenotype.
---

# Disease LoF research (PLoTM for a new disease)

You are running the PLoTM pipeline — *which transcriptional programs of a cell type carry a disease's
rare loss-of-function (LoF) genetic burden, and which regulators drive them* (Ota et al. 2026) — for a
new disease. The heavy pipeline lives in the **`proglof-mcp`** connector; disease background, genetics,
phenotype IDs, and the single-cell substrate come from Claude Science's **Featured connectors**. Your job
is to orchestrate them correctly and interpret the output honestly.

Work inside the granted project folder: write everything under `research/<disease>/` and keep a running
`research/<disease>/brief.md` (copy `project/research_brief_template.md`).

## The loop

1. **Background.** Resolve the disease to an ontology id (OLS → MONDO/EFO), then use literature access +
   OpenAlex/PubMed to write 5–8 sentences on its immune/CD4-T biology and known genetic architecture
   (common vs rare variant). Note the *a priori* programs you'd expect to carry signal (e.g. type-I IFN
   for lupus). Record sources.
2. **Phenotype panel.** Using GWAS Catalog / FinnGen / Open Targets, assemble **3–6 phenotype facets from
   BOTH GeneBass and Backman** where possible (complementary case definitions raise robustness). Each
   needs a source + identifier: a **GeneBass phenocode** or a **Backman GWAS-Catalog accession**. Keep an
   endophenotype (e.g. an autoantibody) as its own trait group if one exists. → `write_config`.
3. **Substrate.** Pick the CRISPRi/CRISPRa Perturb-seq atlas of the relevant cell type (CELLxGENE /GEO),
   set `perturbseq.url`, `pseudobulk_h5ad`, and — critically — `obs_map`. If it has no stimulation axis,
   set one condition. (Reuse the reference CD4-T atlas if the disease is T-cell-mediated.)
4. **Genetics.** `fetch_genebass_burden` + `fetch_backman_burden`, then `run_genebayes`. **Then call
   `qc_report` and read it before trusting anything** — flagged traits are produced but excluded from
   interpretation.
5. **Programs (LoF-independent, the expensive part, once per condition).** `run_phase input` → `cnmf` →
   `stageB`. Routes to the SSH cluster if configured.
6. **Burden + figures + analysis.** `run_phase thresholds` → `burden` → `figures` → `annotate` →
   `analysis`. Then `read_results` (meta, tidy, summary) and write the report.

Resume freely — `run_phase` skips completed work.

## Guardrails — do not violate (see `reference/guardrails.md` for the full statements)

- **γ is condition-independent** — one value per gene per trait; there is no Δγ and no donor-level
  G×condition model. Condition effects live only in the programs/β. Compare conditions at the program
  level (top-gene Jaccard), never at the gene-γ level.
- **GeneBass and Backman are the same UK Biobank cohort** — not independent. Report within- and
  cross-source separately; the meta-analysis is random-effects with I²/Q shown.
- **Never share a fixed |γ| cutoff across sources** (Backman ≈ 60× GeneBass). Use `lof_threshold: q99`.
- **Report the permutation P, not Fisher, for Fig 5a.**
- **Honor the QC gates** (`PRIOR_DOMINATED` / `SIGN_DEGENERATE` / `BORDERLINE`) — flag, don't drop, and
  exclude from robust claims.
- **cNMF is independent per condition** — rest-P2 ≠ stim-P2.
- Pseudobulk Stage-B gate is `stageB_min_profiles: 4`, not the single-cell 10.

If a result contradicts a guardrail (e.g. a "reversal" that is just the one nominally-significant
phenotype), say so plainly rather than reporting the headline.

## Which connector for which step

`reference/connectors.md` maps every step to the exact Featured connector (OLS, MyGene, GWAS Catalog,
FinnGen, Open Targets, literature access, OpenAlex, PubMed, CELLxGENE, GEO) and to the `proglof-mcp`
tools. `reference/pipeline.md` lists the phases, what each produces, and the verify gate after each.

## Deliverable

Write `research/<disease>/report.md` following `reference/report_template.md`: background + sources, the
phenotype panel with QC, the robust cross-source programs (meta FDR, I², LOO, both sources agree), the
phenotype/source-specific ones, the excluded/flagged ones, and an explicit "what the data cannot support"
section. Save figures as artifacts so the reviewer can check every numeric claim against `read_results`.
