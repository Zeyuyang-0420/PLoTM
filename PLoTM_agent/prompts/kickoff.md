# Kickoff prompt

Paste into a Claude Science session in the project (the `disease-lof-research` skill loads automatically).
Replace `<DISEASE>`.

---

Research **<DISEASE>** end-to-end with the PLoTM pipeline.

1. Build the background: resolve <DISEASE> to a MONDO/EFO id, then summarise its CD4⁺ T-cell biology and
   genetic architecture (common vs rare variant) from the literature, and note which programs you expect
   a priori to carry LoF signal. Cite sources.
2. Assemble a LoF phenotype panel of 3–6 facets from **both GeneBass and Backman**, with verified
   phenocodes/accessions from GWAS Catalog / Open Targets. Keep any endophenotype as its own group.
3. Choose the CRISPRi Perturb-seq substrate (reuse the reference CD4-T atlas if the disease is
   T-cell-mediated) and `write_config`.
4. Run the pipeline via the proglof connector: fetch burden (both sources) → run_genebayes → **qc_report
   (read it before trusting anything)** → input → cnmf → stageB → thresholds → burden → figures →
   annotate → analysis. Route heavy phases to the cluster. Resume as needed.
5. Read the meta/tidy/summary results and write `research/<disease>/report.md` following the skill's
   report template: robust cross-source programs vs phenotype-specific vs QC-excluded, plus an explicit
   "what the data cannot support" section.

Respect the guardrails throughout (γ is condition-independent; GeneBass and Backman are one cohort; per-
trait q99 threshold; permutation P not Fisher; honor QC flags). Ask me to approve each new folder, host,
and job. When done, run the reviewer over the report.

---

Follow-ups that also trigger the skill:
- "Which CD4 T-cell programs carry the <DISEASE> LoF burden, and are they robust across sources?"
- "Add the Backman seronegative facet and re-run just the burden + analysis phases."
- "Compare <DISEASE> to the RA reference — which robust programs are shared?"
