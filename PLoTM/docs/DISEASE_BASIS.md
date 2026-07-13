# Optional disease single-cell cNMF basis

By default PLoTM is a **pure Perturb-seq** pipeline (tag `proglof-v1.0-pure-perturbseq`): cNMF gene
programs are learned from the Perturb-seq pseudobulk itself. This document describes the **optional**
mode (v1.1) where programs are instead learned from **disease single-cell data** (e.g. patient CD4⁺
T cells) and the Perturb-seq conditions are *projected* onto that fixed basis.

Turn it on by adding a `disease_singlecell:` block to the config; remove/comment it for pure mode.
Everything downstream (Stage B, Stage C, Fig 4C, Fig 5a) is unchanged — it just consumes different
usages + program identities.

```
                       DEFAULT (v1.0)                         OPTIONAL disease basis (v1.1)
  programs      cNMF on Perturb-seq pseudobulk         cNMF on DISEASE single-cell data (K=30)
  usages        Perturb-seq's own cNMF usages          Perturb-seq PROJECTED onto the fixed disease basis
  meaning       "Perturb programs -> trait"            "Perturb regulators -> DISEASE programs -> trait"
```

## Data flow (disease mode)

```
disease single-cell (.rds/.h5ad, raw counts)
        │  03b/rds_to_counts.R  +  03b/prep_disease_singlecell.py
        │     • intersect genes with the Perturb-measured set   • Stage-A QC (filter_cells)
        ▼
   data/perturbseq/filtered_data/<disease_ct>.h5ad
        │  03_cnmf_stageA/run_cnmf.sh  (K = disease_singlecell.K, e.g. 30)
        ▼
   DISEASE cNMF basis:  gene_spectra_score / gene_spectra_tpm (K × genes) + usages
        │
        │  for EACH Perturb-seq condition:
        │  03b/project_onto_disease.py — X = TPM(condition)[:, disease genes];  U = argmin_{U≥0} ‖X − U·W‖²
        │  (W = disease gene_spectra_tpm, FIXED, update_H=False, assert-checked)
        ▼
   results/cNMF/<CT>/test1/test1.{usages,gene_spectra_score,gene_spectra_tpm}.k_<K>.dt_0_4.*
   results/cNMF/<CT>/{program_estimability.tsv, alive_programs.txt}
        │
        ▼  (unchanged) Stage B β_x(P) → Stage C burden → Fig 4C → Fig 5a
```

## Config

Add to `config/<disease>.yaml` (see `config/RA_singlecell.yaml` for a ready RA example):

```yaml
disease_singlecell:
  rds:   "data/disease_sc/RA_CD4T.rds"   # Seurat .rds (or h5ad: raw counts in .X / a layer)
  counts_layer: counts.1                 # assay layer holding raw counts
  gene_id: symbol                        # symbol | ensg
  K: 30                                  # disease-basis program count (downstream K becomes this)
  intersect_perturb_measured: true       # restrict to Perturb-measured genes (projectable space)
  ct_name: RA_CD4T
```

Run exactly as pure mode — the driver auto-detects the block:

```bash
PY=.../cnmf-env/python CNMF=.../cnmf RSCRIPT=.../Rscript \
  drivers/run_disease.sh config/RA_singlecell.yaml all
# or just the program-basis rebuild + downstream:
  drivers/run_disease.sh config/RA_singlecell.yaml cnmf
  drivers/run_disease.sh config/RA_singlecell.yaml stageB
  drivers/run_disease.sh config/RA_singlecell.yaml burden
  drivers/run_disease.sh config/RA_singlecell.yaml figures
```

## Estimability (the one behaviour disease mode adds)

Stage B contrasts each perturbation against the non-targeting controls. A projected program whose
usage has **zero variance among the NTC controls** admits no contrast, so it has **no β_x(P)**. Such
programs are:

- detected in `project_onto_disease.py` → written to `results/cNMF/<CT>/{program_estimability.tsv,
  alive_programs.txt}`;
- guarded in `regulatory_effectsize.R` (emits an empty β file, no crash);
- **auto-excluded** from the Stage C burden merge and the Fig 5 regsubsets candidate set — the
  scripts auto-discover `alive_programs.txt` per CT (env `ALIVE_PROGRAMS` overrides). In pure
  Perturb-seq mode that file is absent, so all K programs are used and behaviour is **identical**.

`gene_spectra_score` (program identity) is the disease's and shared across conditions; only the
usages differ per condition. So program-burden (Fig 4C x-axis) is condition-independent by
construction — only the regulator axis (β) changes between conditions.

## New / changed files

```
config/disease.example.yaml            + optional disease_singlecell block (commented)
config/RA_singlecell.yaml              ready RA example (disease basis, K=30)
src/03b_disease_basis/
  rds_to_counts.R                      Seurat .rds -> portable counts bundle
  prep_disease_singlecell.py           build cNMF input h5ad (gene intersect + QC)
  project_onto_disease.py              fixed-basis NNLS projection (per condition)
  run_disease_basis.sh                 orchestrates prep -> cNMF -> project
src/07_figures/fig5a_regulator_model.R paper-style regulator->program->trait Fig 5a
drivers/run_disease.sh                 auto-detects disease_singlecell; K=disease K; passes K to Stage B
src/04_regulatory_stageB/regulatory_effectsize.R   + optional K & min-profiles args, NTC guard (default 60 = pure mode)
src/05_burden_stageC/burden_program_regulators.R   + auto-discovered alive-programs restriction (no-op in pure mode)
src/07_figures/fig5_step1_permutation.R            + same alive-programs restriction
```

## Validation status (honest)

The individual stages were validated bit-for-bit against the RA reference already computed under
`tcell_perturbseq/cNMF_RA_analysis/RA_CD4T_external/`:

- **projection** — `project_onto_disease.py` reproduces the reference usages to max |Δ| = 1.2e-15;
  estimability call matches (25/30 for Stim48hr, non-estimable {2,9,14,15,29}).
- **Stage B (K-parameterised + guard)** — reproduces 9,657-gene β files; emits empty for
  non-estimable programs without crashing.
- **Fig 5a regulator model** — edges byte-identical to the reference figure.

A full **from-scratch** run (download the Perturb-seq atlas + GeneBayes GPU posteriors + cNMF on the
disease data + all conditions × traits) is not executed in CI here — it needs the large downloads and
a GPU. The RA disease basis and its Fig 4C / Fig 5a live under `RA_CD4T_external/` as the worked
reference.
