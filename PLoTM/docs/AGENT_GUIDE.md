# Agent guide — running PLoTM for a new disease

This pipeline is designed to be driven by an autonomous agent primarily from a new
`config/<disease>.yaml`. The upstream stages (`input · cnmf · stageB · thresholds`) are fully
config-driven; the reference `burden · figures · analysis` stages run as-is for RA and take new paths
via env vars; only the `download`/`posterior` fetch scripts need a per-disease edit (you must supply
each phenotype's source identifier — burden stats can't be auto-discovered from a name). See the
INTEGRATION STATUS block at the top of `drivers/run_disease.sh`. This guide is the decision procedure,
plus the **non-obvious data facts** that an agent must respect or it will silently produce wrong results.

## The 6-step loop

1. **Pick the Perturb-seq substrate.** Any CRISPRi/CRISPRa Perturb-seq atlas of the relevant cell type
   with a pseudobulk profile matrix (one profile per run×donor×condition×guide). Set
   `perturbseq.url`, `perturbseq.pseudobulk_h5ad`, and — critically — `perturbseq.obs_map` so the
   dataset's obs columns map to `{perturbed_gene, ntc_flag, donor, run, condition, keep_flag}`. If the
   dataset has no stimulation axis, set `conditions:` to the single condition it has.
2. **List the disease phenotypes × source.** Give 3–6 phenotype facets from **both** GeneBass and
   Backman where possible (complementary case definitions raise robustness). Each needs its source and
   identifier (GeneBass phenocode / Backman GWAS-Catalog accession). Keep an endophenotype (e.g. an
   autoantibody) as its own group if the disease has one.
3. **Run `download` → `posterior`.** Fetch burden stats, fit GeneBayes γ. Then **read the QC table**
   (`qc_posteriors.py`) before trusting anything — see the gates below.
4. **Run `input` → `cnmf` → `stageB`** once per condition (LoF-independent, the expensive part).
5. **Run `thresholds` → `burden` → `figures`** — the per-(trait,condition) matrix. Resumable.
6. **Run `annotate` → `analysis`.** Read the meta-analysis and permutation outputs; write conclusions
   that separate robust cross-source findings from phenotype-specific or QC-unstable ones.

## Data facts an agent MUST respect (these are the traps)

1. **γ is condition-independent.** The GeneBayes posterior is human population genetics — one value per
   gene per phenotype, reused for every condition. So **Δγ between conditions ≡ 0**, and there is **no
   donor-level `Y ~ G + condition + G×condition` model to fit** (the genetics and the Perturb-seq share
   no individuals). Condition-dependence lives ONLY in the cNMF programs and β_x(P). Compare conditions
   at the **program level** (matched programs), never at the gene-γ level.
2. **Never share a fixed |γ| threshold across sources.** Backman γ is ~60× the GeneBass scale (sd
   0.17–0.51 vs 0.005–0.008). A fixed cutoff selects ~1% of genes for GeneBass but 60–100% for Backman,
   emptying the Fig-5a background. Use `lof_threshold: q99` (the 99th percentile of |γ| per trait) — for
   a GeneBass trait it lands near the paper's own values; for Backman it scales correctly.
3. **GeneBass and Backman are NOT independent cohorts** (both are UK Biobank exomes). Report
   within-source and cross-source evidence separately; the meta-analysis pools with random effects and
   reports I² / Cochran Q so shared-cohort inflation is visible.
4. **Honor the QC gates — flag, don't drop.** `qc_posteriors.py` flags a trait
   `PRIOR_DOMINATED` if corr(post_mean, prior_mean) > 0.95 (trait data barely move the posterior) and
   `SIGN_DEGENERATE` if <2% or >98% of γ are negative (sign-based tests collapse). Produce these traits
   for completeness but exclude them from meta-analysis and interpretation, and say so.
5. **cNMF programs are condition-specific.** Independent cNMF per condition ⇒ Rest-P*i* ≠ Stim-P*i*.
   Never put them on a shared axis; match across conditions by top-gene Jaccard (reciprocal-best,
   ≥0.10) — the analysis scripts already do this.
6. **Report the permutation P, not the Fisher P, for Fig 5a.** Programs and regulators are selected on
   the same γ vector, so Fisher is strongly anti-conservative. Trust the permutation P.
7. **The Fig-5a regulator step is a 3-of-60 argmax over near-orthogonal predictors.** Individual
   regulator-gene nodes are provisional; lean on program-level meta-analysis and the concordant-gene
   permutation for claims.

## Operational gotchas (will crash or mislead if ignored)

* **`PYTHONNOUSERSITE=1`** always — a broken `~/.local` pandas shadows the conda envs.
* **Run the Fig-5 observed pass before the permutation shards.** Under ~100 concurrent R processes
  `leaps::regsubsets` becomes memory-bandwidth bound (~34 s/perm vs ~9 s) and starves the observed pass
  that gates assembly. `run_batch.sh`/`run_disease.sh` enforce the ordering.
* **Never `pkill -f <pattern>`** when the pattern also matches your own wrapper's command line — it
  kills its own shell (exit 144). Kill by PID selected via `comm==R`.
* **Pseudobulk Stage-B gate is `>4` profiles/gene**, not the paper's single-cell `>10` (pseudobulk has
  ~7 profiles/gene). Set `params.stageB_min_profiles`.
* **A trait with 0 concordant genes** yields a placeholder Fig-5a map (stated, not a crash) — this is
  the honest output for a signless posterior, common for QC-flagged traits.
* **Resumability**: `burden` and `figures` skip completed programs/shards, so a killed run is re-run
  safely at whatever concurrency is free.

## What to verify at each gate

| after | check |
|---|---|
| posterior | QC table: how many traits are flagged? scales sane? |
| cnmf | consensus files exist; k-selection/consensus plots stable |
| stageB | 60 `K60_program*_perturb_effects.txt` per condition |
| thresholds | q99 per trait; flagged traits marked |
| figures | Fig 4C Bonferroni hits; Fig 5a permutation P (not Fisher) |
| analysis | meta I²/LOO stability; overlap permutation (genes vs programs); conclusions separate robust vs unstable |
