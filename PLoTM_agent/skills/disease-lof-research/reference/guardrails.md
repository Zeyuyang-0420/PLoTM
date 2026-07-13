# Guardrails — the traps that silently produce wrong results

These are load-bearing facts about *this* data, not style preferences. Each has burned a real run
(RA and SLE). Violating one yields a confident, wrong answer.

1. **γ is condition-independent.** The GeneBayes posterior is human population genetics — one value per
   gene per phenotype, reused for every condition. Therefore **Δγ between conditions ≡ 0**, and there is
   **no donor-level `Y ~ G + condition + G×condition` model** (the genetics and the Perturb-seq share no
   individuals). All condition-dependence lives in the cNMF programs and their regulatory β. Compare
   conditions at the **program level** (reciprocal-best top-gene Jaccard ≥ 0.10), never at gene-γ.

2. **GeneBass and Backman are NOT independent cohorts** — both are UK Biobank exomes. "Replicates across
   sources" means the same people, analysed twice. Report within-source and cross-source evidence
   separately; the meta-analysis pools with **random effects** and reports **I² / Cochran Q** so
   shared-cohort inflation is visible. Two agreeing sources ≠ two independent replications.

3. **Never share a fixed |γ| threshold across sources.** Backman γ is ~60× the GeneBass scale (sd
   0.17–0.51 vs 0.005–0.008). A fixed cutoff selects ~1% of genes for GeneBass but 60–100% for Backman,
   emptying the Fig-5a background. Always `lof_threshold: q99` (99th percentile of |γ| per trait).

4. **Report the permutation P, not Fisher, for Fig 5a.** Programs and regulators are selected on the same
   γ vector, so Fisher is strongly anti-conservative. The permutation P (NPERM≥5000) is the real test.

5. **Honor the QC gates — flag, don't drop.** `PRIOR_DOMINATED`: corr(post_mean, prior_mean) > 0.95, the
   trait's data barely moved the posterior (low power). `SIGN_DEGENERATE`: <2% or >98% of γ negative, so
   sign-concordance collapses and Fig 5a no longer tests direction. `BORDERLINE`: near a gate. Produce
   these traits for completeness; **exclude them from meta-analysis and interpretation, and say so.**

6. **cNMF programs are condition-specific.** Independent cNMF per condition ⇒ rest-P*i* ≠ stim-P*i*. Never
   put them on a shared axis or claim same-numbered programs are the same program.

7. **The Fig-5a regulator step is a 3-of-60 argmax over near-orthogonal predictors.** Individual
   regulator-gene nodes are provisional. Lean on program-level meta-analysis and the concordant-gene
   permutation for claims. Shared *genes* across phenotypes can exceed chance; shared *programs* often do
   not — check the permutation, don't assume.

8. **A "sign reversal" is usually just the one nominally-significant phenotype.** Before reporting a
   reversal, check all phenotypes' signs and I²: if 4/5 agree and one is merely the only P<0.05, it is not
   a reversal. (This is the RA P31 lesson.)

## Operational

- `PYTHONNOUSERSITE=1` always (a broken `~/.local` pandas shadows the envs).
- Pseudobulk Stage-B gate is `stageB_min_profiles: 4`, not the single-cell 10 (~7 profiles/gene).
- Run the Fig-5 observed pass before the permutation shards; the driver enforces this.
- A trait with 0 concordant genes yields a placeholder Fig-5a map (stated, not a crash) — the honest
  output for a signless posterior, common for QC-flagged traits.
