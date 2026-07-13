# Pipeline methodology — stage by stage

Every stage reads its knobs from `config/<disease>.yaml` and is invoked by
`drivers/run_disease.sh <config> <phase>`. Paths below are relative to the repo root; `<CT>` is a
condition-specific cell-type tag like `RA_rest_pseudobulk`.

## 0 · Download (`src/00_download/`)

* `download_perturbseq.sh` — pseudobulk h5ad for the atlas (one profile per run×donor×condition×guide).
* `fetch_genebass_lof.py`, `fetch_backman_lof.py` — gene-based exome **burden** summary statistics
  (`ensg, gene, chrom, beta, standard_error, z`) for each phenotype, from the two sources, written in a
  common schema. `prepare_lof_inputs_*.py` shape them into GeneBayes response tables.

## 1 · LoF posterior (`src/01_lof_posterior/`, GPU)

`run_genebayes_converged.sh` fits GeneBayes per trait (`--total_iterations 1000`, NGBoost early
stopping ≈ iter 140) → `data/lof/GeneBayes_posterior/<trait>.per_gene_estimates.tsv` with
`post_mean` (γ), `prior_mean`, and 95% CI. `qc_posteriors.py` computes the **QC gates**
(`PRIOR_DOMINATED`, `SIGN_DEGENERATE`). `ensemble_gamma.py` averages independent converged runs
(the fit is stochastic; r≈0.70 run-to-run; ensembling raises reliability per Spearman–Brown).

## 2 · Input prep (`src/02_input_prep/`)

`make_pseudobulk_input.py` subsets to `keep_flag & condition`, builds the pipeline obs schema
(`gene, gem_group, n_genes, mitopercent, leiden`) from `obs_map`, and writes the raw-count cNMF input
h5ad + covariate metadata. One `<CT>` per condition.

## 3 · Stage A — cNMF programs (`src/03_cnmf_stageA/`)

`filter_cells.py` (mito<0.3, gene/count QC) → `run_cnmf.sh` (prepare/factorize/consensus, K=60, 100
iters, seed 14, 2000 HVGs, dt=0.4). Outputs `gene_spectra_score.k_60.dt_0_4.txt` (60×genes) and
`usages...consensus.txt` (profiles×60). **LoF-independent.**

## 4 · Stage B — regulatory effects (`src/04_regulatory_stageB/`)

`regulatory_effectsize.R` — for each knockdown gene with `>stageB_min_profiles` profiles,
`lm(scale(usage) ~ GROUP + gem_group + n_genes + mitopercent)` vs NTC → β_x(P), p. 60
`K60_program*_perturb_effects.txt` per condition. **LoF-independent.**

## 5 · Stage C — burden (`src/05_burden_stageC/`)

`burden_program_regulators.R` (one program per invocation, parallel). Per program:
**regulator burden** = Pearson r of β_x(P) vs γ + s_het-adjusted `lm(post_mean ~ β + shet)`;
**program burden** = mean γ of the top-100 program genes vs an s_het-stratified permutation null.
Concatenated into `programs_/regulators_enrichment_K60_<trait>.tsv` — the Figure-4C axes.

## 6 · Annotation (`src/06_annotation/`)

`annotate_programs.R` — top-200 genes → `enrichGO(BP)` + hallmark `enricher` + a CD4 marker panel;
universe = expressed genes. `curated_labels.tsv` holds human-checked labels used by the networks.

## 7 · Figures (`src/07_figures/`)

* `fig4_burden.R` → **Figure 4C**: x = signed −log10 P of program burden, y = signed −log10 P of
  regulator burden; Bonferroni = 0.05/K, with an adaptive "none significant" note.
* `fig5_step1_permutation.R` (modes `obs` | `perm`) → selected programs/regulators, concordant genes,
  and the permutation null. `fig5_step2_assemble.R` → Fisher P, **permutation P** (report this), and the
  **Figure 5a** regulator→program→trait map (line width = effect strength).
* `multi_trait_network.R` → combined maps across a trait group; conservation-coloured; `ANNOTATED_ONLY=1`
  keeps only labelled programs. `summarize_programs.py` → ranked readout.

## 8 · Condition analysis (`src/08_condition_analysis/`)

Full unthresholded, harmonised comparison of Rest vs Stim (or any condition pair):

* `build_tidy.py` → `tidy_results.tsv` (all programs × traits × conditions, BH-FDR).
* `meta_analysis.py` → DerSimonian–Laird random-effects meta of the regulator β across the disease
  phenotypes (the only program-level effect with an SE); Cochran Q, I², leave-one-phenotype-out,
  within/cross-source.
* `condition_interaction.py` → Jaccard-match programs across conditions, Δβ = β_Stim − β_Rest with
  SE = √(SE²+SE²), and category (constitutive / amplified / attenuated / revealed / sign-reversing).
* `overlap_permutation.py` → does cross-phenotype gene/program overlap exceed chance (per-trait top-hit
  universe null).
* `priority_checks.py`, `make_figures.py` → targeted checks (sign reversals, sub-threshold consistency,
  QC-program correlations) and the heatmap / scatter / forest / network / QC figures.

## Drivers (`drivers/`)

* `run_disease.sh` — config-driven top-level orchestrator (this document's phases).
* `run_batch.sh` — the (trait × condition) burden/figure matrix, resumable, throttled.
* `run_trait.sh` — one (condition, trait) end-to-end; assembles a symlinked `data/` tree so every trait
  reuses the identical cNMF matrix.
* `compute_thresholds.py` — per-trait q99 |γ| threshold + QC flag table.
