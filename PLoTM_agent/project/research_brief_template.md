# Research brief — <DISEASE>

Copy to `research/<disease>/brief.md` and fill as you go. This is the working record the skill maintains
and the reviewer checks against.

## 1. Disease & ontology
- name / synonyms:
- ontology id (MONDO/EFO, via OLS):
- cell type of interest (why): CD4⁺ T cells because …

## 2. Background (literature access + OpenAlex/PubMed)
- immune / CD4-T biology (3–5 sentences):
- genetic architecture — common vs rare variant (2–3 sentences):
- **a-priori programs** expected to carry signal:
- sources (DOIs / accessions):

## 3. Phenotype panel (GWAS Catalog / FinnGen / Open Targets)
| trait name | source | phenocode / accession | class | rationale |
|---|---|---|---|---|
| Genebass_<...> | genebass | <phecode/ICD> | binary/continuous | |
| Backman_2021_<...> | backman | GCST<...> | binary | |
| … | | | | |
- endophenotype (its own group), if any:

## 4. Perturb-seq substrate (CELLxGENE / GEO)
- atlas + accession:
- url / pseudobulk_h5ad:
- obs_map: perturbed_gene= · ntc_flag= · donor= · run= · condition= · keep_flag=
- conditions:

## 5. Config
- `write_config` → path:
- K (60 default; size down if few pseudobulk samples):

## 6. Run log (proglof-mcp)
- posterior + `qc_report`: flagged traits =
- input/cnmf/stageB: done?
- thresholds/burden/figures/annotate/analysis: done?

## 7. Results (`read_results`)
- robust cross-source programs (meta FDR<0.05, I²≈0, LOO-stable, both sources):
- phenotype/source-specific:
- excluded / QC-flagged:
- condition dependence:

## 8. Report
- `research/<disease>/report.md` written? reviewer passed?
