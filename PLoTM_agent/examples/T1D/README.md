# Example — Type 1 diabetes (T1D)

A worked reference for the `disease-lof-research` skill on a disease **distinct from the existing RA and
SLE runs**. It shows the shape of a correct autonomous run and serves as an acceptance target: given the
kickoff prompt with `<DISEASE>=type 1 diabetes`, Claude Science should discover a comparable phenotype
panel, produce a config like `config.yaml`, run the pipeline, and reach the interpretations in `brief.md`.

T1D is chosen because it is **T-cell-mediated** (so the reference CD4⁺ T Perturb-seq atlas is
appropriate) yet has a **different genetic architecture from RA** (HLA-dominated, strong common-variant
component) — a good test of whether the rare-LoF burden side is honestly reported as likely
underpowered, as it was for SLE's IFN axis.

## Files
- `brief.md` — the filled research brief (background, phenotype panel, substrate, expected findings).
- `config.yaml` — the PLoTM config the agent should generate (matches `write_config` output; round-
  trips through `src/cfg_get.py`).

## Acceptance (what a correct run demonstrates)
1. A phenotype panel with **both** GeneBass and Backman facets and verified identifiers.
2. `qc_report` read and flagged traits excluded from the robust section.
3. A report that **separates** any robust cross-source program from phenotype-specific ones, reports the
   **permutation** P for Fig 5a, and includes an explicit "what the data cannot support" section.
4. Honest handling if the rare-LoF side is underpowered (state it; don't manufacture a hit) — the same
   discipline the RA and SLE runs used.

> Identifiers below marked *(verify)* are canonical phecodes/ICD codes; Backman GWAS-Catalog accessions
> must be confirmed in the GWAS Catalog / Open Targets connector at run time (do not hard-code an
> unverified accession).
