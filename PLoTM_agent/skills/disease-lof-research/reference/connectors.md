# Step → connector map

Claude Science's **Featured connectors** do discovery; the **proglof-mcp** connector does the pipeline.
Name the source in your request or let Claude pick; queries appear as expandable steps.

| step | use | connector(s) |
|---|---|---|
| disease → ontology id | resolve name to MONDO/EFO; get associated genes | **OLS** (Genes & Ontologies), **MyGene** |
| immune / CD4-T biology background | full-text of key reviews & primary papers | **literature access** (built-in) + **OpenAlex/arXiv** (Literature Graph), **PubMed** (directory), **bioRxiv** |
| genetic architecture (common vs rare) | which loci/genes, effect classes, evidence | **GWAS Catalog**, **FinnGen**, **BioBank Japan** (Human Genetics); **Open Targets** (Clinical Genomics) |
| candidate LoF phenotypes + identifiers | find GeneBass phenocodes & Backman accessions for the disease's facets | **GWAS Catalog** (Backman accessions), **Open Targets** (trait mapping); GeneBass phenocodes from the disease's phecodes/ICD |
| variant / gene sanity checks | is a top-γ gene a known LoF-intolerant / disease gene? | **gnomAD**, **ClinVar** (Variants); **Open Targets**, **ClinGen** |
| Perturb-seq / patient single-cell substrate | choose the cell-type atlas; find its accession | **CELLxGENE / CellGuide**, **GEO**, **ArrayExpress** |
| program annotation cross-check | what is this program's biology | **GO / Reactome** (Genes & Ontologies), **Human Protein Atlas**, **STRING** |

## proglof-mcp tools (the pipeline)

`write_config` → `fetch_genebass_burden` / `fetch_backman_burden` → `run_genebayes` → `qc_report` →
`run_phase input|cnmf|stageB|thresholds|burden|figures|annotate|analysis` → `read_results`.

Set the two `fetch_*` tools and `run_genebayes` to *Ask each time*; `write_config`, `read_results`,
`qc_report`, `run_phase` are safe to *Always allow*.

## What proglof-mcp does NOT do

It does not discover phenotype identifiers, fetch papers, or map ontologies — those are the Featured
connectors above. `write_config` only records the identifiers you discovered; garbage in → garbage out,
so verify each phenocode/accession against GWAS Catalog / Open Targets before writing it.
