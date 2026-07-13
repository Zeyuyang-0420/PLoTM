# Report template — `research/<disease>/report.md`

Write for a skeptical immunogeneticist. Every number must be backed by a `read_results` call the reviewer
can re-check. Separate robust from unstable; never bury a caveat.

```markdown
# <Disease>: CD4⁺ T-cell programs vs rare-LoF burden

## Background
<5–8 sentences: immune/CD4-T biology, genetic architecture (common vs rare variant), and the programs
you expected a priori to carry signal. Cite sources (DOIs / accessions).>

## Data
- Perturb-seq substrate: <atlas, accession, conditions>
- LoF phenotypes (panel): table of {trait, source, phenocode/accession, class, QC flag}
- Sources are BOTH UK Biobank exomes (GeneBass, Backman) — not independent cohorts.

## Robust, cross-source programs  (the trustworthy result)
Table from `read_results(kind="meta")`: program | annotation | condition | pooled β | FDR | I² |
cross-source | LOO-stable. Keep rows with FDR<0.05, I²≈0, LOO-stable, both sources agreeing. State the
biology (does it match the a-priori expectation?).

## Phenotype- or source-specific signals  (weaker)
Programs significant in one source/phenotype only, or with high I². Label them as such.

## Condition dependence
Program-level (matched by top-gene Jaccard). Expect mostly constitutive because γ is
condition-independent. Note any genuine rest/stim-specific program.

## Excluded / QC-flagged
List PRIOR_DOMINATED / SIGN_DEGENERATE / BORDERLINE traits and the noise/QC programs, with the reason
each is excluded (from `qc_report`).

## What the data cannot support
Explicitly: no Δγ, no donor G×condition, no independent replication (single cohort), Fig-5a regulator
nodes provisional. Say what would be needed to go further.

## Figures
Fig 4C, Fig 5a (permutation P), meta forest, condition heatmap — saved as artifacts.
```

## Reviewer checklist (what Claude Science's reviewer should be able to confirm)

- Every meta β/FDR/I² in the report matches `read_results(kind="meta")`.
- Every "significant" uses the **permutation** P (Fig 5a) or FDR (meta), never Fisher.
- No QC-flagged trait appears in the robust section.
- No cross-condition claim treats same-numbered programs as identical.
- The phenotype identifiers in the config match what GWAS Catalog / Open Targets returned.
