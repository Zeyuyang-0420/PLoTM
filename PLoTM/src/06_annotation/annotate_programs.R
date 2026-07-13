#!/usr/bin/env Rscript
# annotate_programs.R <CT> [TOPN=200] [OUTDIR=program_annotation]
#
# Annotates each of the 60 cNMF programs of one condition, following the paper's Methods
# ("Annotation of programs to biological pathways"), which uses three orthogonal lines of evidence:
#
#   (1) GO / MSigDB-hallmark enrichment of the TOP-200 loading genes, using clusterProfiler's
#       enrichGO() and enricher(), with EXPRESSED GENES as the background/universe.   <-- implemented
#   (2) Transcription-factor binding-site enrichment from ENCODE ChIP-seq in K562.    <-- NOT
#       transferable: those ChIP scores are K562-specific and the paper's 167-TF panel was
#       validated against K562 knockdown DEGs. Omitted, and flagged as such.
#   (3) Co-expression with marker genes for predefined cell types / pathways.         <-- implemented
#       as marker-set overlap of the top-200 genes, using a CD4+ T-cell marker panel in place of
#       the paper's erythroid/myeloid/K562 panels.
#
# The paper's curation rule: if the program's top-loading genes contain cell-type marker genes,
# annotate it as that cell type/state; otherwise take the most significant non-ambiguous GO/hallmark
# term from the top 10. We emit the evidence; the final label is curated in curate_annotations.py.
#
# Genes are ranked by gene_spectra_score (the same ranking the burden test and the Fig-5 map use to
# define "program genes"). NB the Methods text says "non-negative loadings"; gene_spectra_tpm is also
# written for comparison.

suppressMessages({library(clusterProfiler); library(org.Hs.eg.db)})
args <- commandArgs(trailingOnly = TRUE)
CT     <- args[1]
TOPN   <- if (length(args) >= 2) as.numeric(args[2]) else 200
OUTDIR <- if (length(args) >= 3) args[3] else "program_annotation"
K <- 60
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

RESULTS <- Sys.getenv("RESULTS_ROOT", "/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis")
DATA    <- Sys.getenv("DATA_ROOT", "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data")

GEP <- read.table(file.path(RESULTS, "cNMF", CT, "test1",
                            paste0("test1.gene_spectra_score.k_", K, ".dt_0_4.txt")),
                  header = TRUE, stringsAsFactors = FALSE)
GEP <- t(GEP); colnames(GEP) <- paste0("P", 1:K)
corresp <- read.table(file.path(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID"), header = FALSE,
                      stringsAsFactors = FALSE)
corresp <- corresp[!duplicated(corresp[, 1]), ]; row.names(corresp) <- corresp[, 1]
sym <- corresp[row.names(GEP), 2]
keep <- !is.na(sym) & !duplicated(sym)
GEP <- GEP[keep, ]; row.names(GEP) <- sym[keep]
universe <- row.names(GEP)                      # expressed-gene background (paper's approach)
cat(sprintf("[annot] %s: %d programs x %d expressed genes (universe)\n", CT, ncol(GEP), length(universe)))

## hallmark gene sets (MSigDB) if msigdbr available
H <- NULL
if (requireNamespace("msigdbr", quietly = TRUE)) {
  hh <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
  if (is.null(hh)) hh <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", collection = "H"), error = function(e) NULL)
  if (!is.null(hh)) {
    gs <- if ("gene_symbol" %in% names(hh)) hh$gene_symbol else hh$gs_symbol
    nm <- if ("gs_name" %in% names(hh)) hh$gs_name else hh$gs_id
    H <- data.frame(term = nm, gene = gs, stringsAsFactors = FALSE)
    cat(sprintf("[annot] hallmark sets: %d terms\n", length(unique(H$term))))
  }
}
if (is.null(H)) cat("[annot] msigdbr unavailable -> GO only\n")

## CD4+ T-cell marker panel (paper's line 3, retargeted from K562/erythroid to T cells)
markers <- list(
  Naive_CM        = c("CCR7","SELL","TCF7","LEF1","IL7R","MAL","NOSIP"),
  Th1             = c("TBX21","IFNG","CXCR3","IL12RB2","IL18R1","CCL5"),
  Th2             = c("GATA3","IL4","IL5","IL13","IL17RB","PTGDR2","HPGDS"),
  Th17            = c("RORC","CCR6","IL17A","IL17F","IL23R","IL1R1","KLRB1"),
  Treg            = c("FOXP3","IL2RA","IKZF2","CTLA4","TNFRSF18"),
  Tfh_Tph         = c("PDCD1","CXCR5","CXCL13","ICOS","BTLA","MAF"),
  Cytotoxic       = c("GZMK","GZMB","GZMH","GZMA","NKG7","PRF1","EOMES","CST7","CTSW","KLRG1","KLRD1"),
  Exhaustion      = c("PDCD1","LAG3","TIGIT","HAVCR2","CD160","TOX"),
  CellCycle_S     = c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","TYMS","RRM2","GINS2"),
  CellCycle_G2M   = c("CDK1","CCNB1","CCNB2","TOP2A","MKI67","UBE2C","BUB1","PLK1","AURKB"),
  Interferon      = c("ISG15","MX1","MX2","OAS1","OAS2","OAS3","IFIT1","IFIT2","IFIT3","STAT1","IRF7","HERC6"),
  Translation_Ribo= c("RPL3","RPL4","RPS3","RPS6","EIF4G1","EIF3L","NCL","NOP16","NPM1","TCP1","CCT6A"),
  OXPHOS_Mito     = c("NDUFA13","NDUFB8","COX5A","ATP5F1A","MRPL16","MRPL42","UQCRQ","SDHB","ATP5PO"),
  Activation_Costim=c("TNFRSF4","TNFRSF9","TNFRSF18","CD69","IL2RA","ICOS","CD40LG"),
  Apoptosis_p53   = c("BAX","CDKN1A","MDM2","ZMAT3","FDXR","BCL2L11","TP53I3","PHLDA3","TRIM22"),
  Epithelial      = c("EPCAM","CDH1","KRT8","KRT18","DSG2","ESRP1","SPINT1","SPINT2","LSR","AP1M2"),
  Cytokine_Inflam = c("IL6","IL1B","IL1R1","TNF","CXCL8","PTGS2","CCL20"),
  Prostaglandin   = c("PTGS2","PTGDS","HPGDS","PLA2G4A","PTGDR2","PTGER2"),
  RANKL_Bone      = c("TNFSF11","TNFRSF11A"),
  TCR_identity    = c("CD4","CD3E","CD3D","CD247","LCK","ZAP70","THEMIS","LAT","ITK","IL2RG")
)

rows <- list()
for (i in 1:K) {
  p <- paste0("P", i)
  top <- names(sort(GEP[, p], decreasing = TRUE))[1:TOPN]

  ego <- tryCatch(enrichGO(gene = top, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                  error = function(e) NULL)
  gdf <- if (!is.null(ego)) as.data.frame(ego) else data.frame()
  go_terms <- if (nrow(gdf)) paste(sprintf("%s (q=%.1e)", gdf$Description[1:min(5, nrow(gdf))],
                                           gdf$p.adjust[1:min(5, nrow(gdf))]), collapse = " | ") else "-"

  hall <- "-"
  if (!is.null(H)) {
    eh <- tryCatch(enricher(gene = top, universe = universe, TERM2GENE = H,
                            pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
    hdf <- if (!is.null(eh)) as.data.frame(eh) else data.frame()
    if (nrow(hdf)) hall <- paste(sprintf("%s (q=%.1e)", hdf$Description[1:min(3, nrow(hdf))],
                                         hdf$p.adjust[1:min(3, nrow(hdf))]), collapse = " | ")
  }

  mk <- sapply(markers, function(g) sum(g %in% top))
  mk <- mk[mk > 0]
  mk_s <- if (length(mk)) paste(sprintf("%s:%d/%d", names(mk), mk, sapply(markers[names(mk)], length)),
                                collapse = ",") else "-"

  rows[[i]] <- data.frame(CT = CT, Program = p,
                          n_GO = nrow(gdf),
                          top_GO = go_terms,
                          top_hallmark = hall,
                          marker_hits = mk_s,
                          top25_genes = paste(top[1:25], collapse = ","),
                          stringsAsFactors = FALSE)
  cat(sprintf("  %s: %d GO terms | markers: %s\n", p, nrow(gdf), substr(mk_s, 1, 60)))
}
out <- do.call(rbind, rows)
f <- file.path(OUTDIR, paste0(CT, "_program_annotation_top", TOPN, ".tsv"))
write.table(out, f, row.names = FALSE, sep = "\t", quote = TRUE)
cat(sprintf("[annot] DONE -> %s\n", f))
