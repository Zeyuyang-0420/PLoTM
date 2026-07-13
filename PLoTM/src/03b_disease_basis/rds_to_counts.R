#!/usr/bin/env Rscript
# rds_to_counts.R <rds> <counts_layer> <out_prefix>
#
# Export a Seurat .rds disease single-cell object to a portable counts bundle that
# prep_disease_singlecell.py can assemble into an AnnData without needing R in the python env:
#   <out_prefix>.counts.mtx    genes x cells, raw integer counts (MatrixMarket)
#   <out_prefix>.genes.txt     gene identifiers (rownames)
#   <out_prefix>.barcodes.txt  cell barcodes (colnames)
#   <out_prefix>.metadata.csv  obs (cell_barcode + all meta.data columns; percent.mt kept if present)
#
# counts_layer names the assay layer holding RAW counts (Seurat v5: e.g. "counts" or "counts.1").

suppressMessages({library(Seurat); library(SeuratObject); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
rds <- a[1]; layer <- if (length(a) >= 2 && nzchar(a[2])) a[2] else "counts"
out <- a[3]
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(rds)
stopifnot(inherits(obj, "Seurat"))
da <- DefaultAssay(obj)
m <- tryCatch(LayerData(obj, assay = da, layer = layer),
              error = function(e) SeuratObject::GetAssayData(obj, assay = da, slot = "counts"))
m <- as(m, "CsparseMatrix")
cat(sprintf("[rds] %s assay=%s layer=%s : %d genes x %d cells\n", basename(rds), da, layer, nrow(m), ncol(m)))
stopifnot(all(m@x == round(m@x)))                     # raw integer counts

Matrix::writeMM(m, paste0(out, ".counts.mtx"))
writeLines(rownames(m), paste0(out, ".genes.txt"))
writeLines(colnames(m), paste0(out, ".barcodes.txt"))
md <- obj@meta.data; md$cell_barcode <- rownames(md)
md <- md[, c("cell_barcode", setdiff(colnames(md), "cell_barcode"))]
write.csv(md, paste0(out, ".metadata.csv"), row.names = FALSE)
cat(sprintf("[rds] wrote %s.{counts.mtx,genes.txt,barcodes.txt,metadata.csv}\n", out))
