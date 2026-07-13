# Research brief — Type 1 diabetes (T1D)

## 1. Disease & ontology
- name / synonyms: type 1 diabetes mellitus, T1D, insulin-dependent diabetes
- ontology id: **MONDO:0005147** / EFO:0001359 (via OLS)
- cell type: CD4⁺ T cells — T1D is a T-cell-mediated autoimmune destruction of pancreatic β-cells;
  autoreactive CD4⁺ T help and Treg dysfunction are central.

## 2. Background
- CD4-T biology: islet-antigen-specific CD4⁺ T cells drive β-cell autoimmunity; Treg
  insufficiency and an IL-2/STAT5 signalling defect are recurrent themes; type-I/II IFN and cytotoxic
  CD4 programs are implicated in progression.
- genetic architecture: **HLA class II (DR/DQ) dominates**; strong common-variant component (INS, PTPN22,
  IL2RA, CTLA4, IFIH1, …). The **rare coding-LoF** contribution is comparatively small — so, as for SLE,
  expect the rare-LoF burden side to be underpowered and possibly null; report honestly.
- a-priori programs: Treg / IL-2-STAT5, cytotoxic CD4 (GZMK/EOMES), IFN response, TCR/CD4 identity.
- sources: MONDO/EFO (OLS); Open Targets T1D (MONDO:0005147); reviews via literature access + PubMed.

## 3. Phenotype panel
| trait name | source | phenocode / accession | class | rationale |
|---|---|---|---|---|
| Genebass_T1D          | genebass | ICD-10 **E10** *(verify)*        | binary | primary T1D case definition |
| Genebass_T1D_phecode  | genebass | phecode **250.1** *(verify)*    | binary | complementary phecode definition |
| Genebass_T1D_early    | genebass | early-onset / IDDM subset       | binary | onset-stratified facet |
| Backman_2021_T1D      | backman  | GWAS-Catalog **GCST……** *(verify)* | binary | independent case definition, other source |
| Backman_2021_T1D_alt  | backman  | GWAS-Catalog **GCST……** *(verify)* | binary | alternate Backman T1D definition |
- endophenotype group: (optional) GAD/IA-2 autoantibody trait if a burden facet exists — its own group.

## 4. Perturb-seq substrate
- atlas: reuse the reference genome-scale **CD4⁺ T-cell CRISPRi Perturb-seq (GWCD4i)** — T1D is
  T-cell-mediated, so the same substrate + Stage-B β apply.
- conditions: Rest, Stim48hr (reuse the reference cNMF programs; only the γ vectors change).
- obs_map: as the RA reference.

## 5. Config
- `write_config` → `config/T1D.yaml` (see `config.yaml`), K=60, `lof_threshold: q99`,
  `stageB_min_profiles: 4`.

## 6. Expected run notes
- `qc_report`: watch for PRIOR_DOMINATED among the lower-power Backman facets (small case counts) — as in
  RA's seropos/seroneg facets.
- programs/Stage-B are reused from the reference (LoF-independent) — the expensive `cnmf`/`stageB` need
  not be re-run if pointing at the reference CD4-T matrix; only `thresholds`→`analysis` per new γ.

## 7. Expected findings (acceptance)
- **Likely outcome:** no robust cross-source rare-LoF program (HLA/common-variant disease; rare-LoF
  underpowered). Any Bonferroni pass on a *borderline* housekeeping program is an artefact — ring it and
  exclude, exactly as SLE's P14 and RA's Rest P3.
- If a program does survive (meta FDR<0.05, I²≈0, LOO-stable, both sources), report it with its biology
  (plausible candidates: Treg / IL-2-STAT5, cytotoxic CD4).
- **What the data cannot support:** no Δγ / donor G×condition; single cohort (no independent replication);
  Fig-5a regulator nodes provisional.

## 8. Report
- `research/type-1-diabetes/report.md` per the skill template; reviewer confirms every meta number.
