#!/usr/bin/env Rscript
# plot_multi_trait_network.R <CT> <GROUP_NAME> <TRAIT1,TRAIT2,...> [FIGDIR] [LABELS] [NPERM_SIGN]
#
# Generalisation of plot_two_trait_network.R to N traits, to surface the CONSERVED programs and
# KEY (shared) regulator genes across a group of related LoF traits in one condition.
#
#   left   : union of concordant top-|gamma| genes across the group's traits
#            (each trait uses ITS OWN |gamma| threshold = 99th pct, the same rule as its Fig5a map),
#            coloured by CONSERVATION = the number of traits in which the gene is concordant.
#   middle : programs selected (as a program OR a regulator) by any trait, labelled with the curated
#            annotation, coloured by how many traits selected them.
#   right  : one node per trait (short labels), spread across the gene axis.
#
# LINE WIDTH == STRENGTH, relative WITHIN each edge class (same convention as the two-trait map):
#   gene -> program "regulates" : |beta_x(P)|           (BH<0.05 regulatory effect)
#   gene -> program "member"    : cNMF gene_spectra_score loading (compressed, background)
#   program -> trait            : -log10(P) of the burden statistic that placed it
#                                   regulator programs -> P of the coef in lm(gamma ~ beta_regs + shet)
#                                   program-only        -> shet-matched permutation P of mean gamma (top-TOP)
#
# Per-trait LOF thresholds are read from lof_thresholds.tsv (col `trait`,`thresh`) unless a trait's
# tag file is found directly.  A trait contributing 0 concordant genes (e.g. a degenerate posterior)
# still contributes its program->trait edges; it simply adds no genes to the union.

suppressMessages(library(ggplot2))
args   <- commandArgs(trailingOnly = TRUE)
CT     <- args[1]
GROUP  <- args[2]
TRAITS <- strsplit(args[3], ",")[[1]]
FIGDIR <- if (length(args) >= 4) args[4] else "figures_multi"
LABF   <- if (length(args) >= 5) args[5] else
          "/mnt/scratch/ZY2/Hackathon/RA_cNMF_LoF_pipeline/program_annotation/curated_labels.tsv"
NPERM_SIGN <- if (length(args) >= 6) as.numeric(args[6]) else 2000
F5   <- Sys.getenv("F5_ROOT", "fig5_out")
THF  <- Sys.getenv("THRESH_TSV", "/mnt/scratch/ZY2/Hackathon/pipeline_results/lof_thresholds.tsv")
# ANNOTATED_ONLY=1 : keep only program nodes that carry a curated annotation (declutters the map);
# genes that then have no edge to any shown program are dropped.
ANNOTATED_ONLY <- Sys.getenv("ANNOTATED_ONLY", "0") == "1"
TOP  <- 200; K <- 60
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)
LABEL   <- if (grepl("rest", CT)) "GWCD4i Rest" else "GWCD4i Stim48hr"
RESULTS <- Sys.getenv("RESULTS_ROOT", "/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis")
DATA    <- Sys.getenv("DATA_ROOT", "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data")

## short trait labels (user's naming)
short_of <- function(t) {
  m <- c(Backman_2021_RA_conv = "BK_curated", Backman_2021_RA_M06_conv = "BK_M06",
         Backman_2021_RA_alt_conv = "BK_alt", Genebass_RA_custom_conv = "GB_RA_custom",
         Genebass_RA_M06_conv = "GB_RA_M06", Genebass_RF_conv = "GB_RF")
  if (t %in% names(m)) m[[t]] else sub("_conv$", "", t)
}
## per-trait tag from its LOF threshold
th <- read.table(THF, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thresh_of <- function(t) th$thresh[th$trait == t][1]
tag_of    <- function(t) sprintf("program5_regulator3_LOF%s", format(thresh_of(t), trim = TRUE))

lab <- read.table(LABF, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
lab <- lab[lab$CT == CT, ]; labmap <- setNames(lab$Label, lab$Program)

## ---------------- shared inputs (LoF-independent) ----------------
GEPe <- read.table(file.path(RESULTS, "cNMF", CT, "test1", "test1.gene_spectra_score.k_60.dt_0_4.txt"),
                   header = TRUE, stringsAsFactors = FALSE)
GEPe <- t(GEPe); colnames(GEPe) <- paste0("P", 1:K)
cp <- read.table(file.path(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID"), header = FALSE, stringsAsFactors = FALSE)
cp <- cp[!duplicated(cp[, 1]), ]; row.names(cp) <- cp[, 1]
symv <- cp[row.names(GEPe), 2]; ok <- !is.na(symv) & !duplicated(symv)
GEPs <- GEPe[ok, ]; row.names(GEPs) <- symv[ok]

Breg <- NULL; Preg <- NULL
for (i in 1:K) {
  x <- read.table(file.path(RESULTS, "cNMF_regulation", CT, paste0("K", K, "_program", i, "_perturb_effects.txt")), header = TRUE)
  b <- x[, c(1, 2)]; p <- x[, c(1, 3)]
  colnames(b) <- c("GENE", paste0("P", i)); colnames(p) <- c("GENE", paste0("P", i))
  if (i == 1) { Breg <- b; Preg <- p } else { Breg <- merge(Breg, b, by = "GENE"); Preg <- merge(Preg, p, by = "GENE") }
}
shet <- read.table(file.path(DATA, "shet_10bins.txt"), header = TRUE, stringsAsFactors = FALSE)

## ---------------- per-trait Fig5a artefacts ----------------
edges <- NULL; conc <- list(); selP <- list(); selR <- list(); LOF <- list()
for (t in TRAITS) {
  TAG <- tag_of(t)
  D   <- file.path(t, CT, F5, CT, "P200")                       # per-trait work tree
  if (!dir.exists(D)) D <- file.path(F5, CT, "P200")            # fallback: shared tree
  ef <- file.path(D, paste0(t, "_", TAG, "_map_edges.tsv"))
  # a trait with 0 concordant genes writes a header-only edges file (or none): read.table(header=TRUE)
  # errors on "no lines available", so require >1 line before parsing.
  if (file.exists(ef) && length(readLines(ef, n = 2)) > 1) {
    e <- read.table(ef, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (nrow(e)) edges <- rbind(edges, data.frame(e, trait = t, stringsAsFactors = FALSE))
  }
  cf <- file.path(D, paste0("ConcordantGenes_", t, "_", TAG, ".txt"))
  cg <- tryCatch(read.table(cf, header = TRUE, stringsAsFactors = FALSE)$GENE, error = function(e) character(0))
  conc[[t]] <- if (is.null(cg)) character(0) else cg[!is.na(cg)]
  s <- read.table(file.path(D, paste0("SelectedPrograms_", t, "_", TAG, ".txt")), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  selP[[t]] <- strsplit(s$Program_selected[1], ",")[[1]]; selR[[t]] <- strsplit(s$Regulator_selected[1], ",")[[1]]
  g <- read.table(file.path(DATA, "LoF", "GeneBayes_posterior", paste0(t, ".per_gene_estimates.tsv")),
                  sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
  g <- data.frame(GENE = cp[g$ensg, 2], g); g <- g[!is.na(g$GENE), ]
  g$post_mean[is.infinite(g$post_mean) & g$post_mean > 0] <- max(g$post_mean[!is.infinite(g$post_mean)])
  g$post_mean[is.infinite(g$post_mean) & g$post_mean < 0] <- min(g$post_mean[!is.infinite(g$post_mean)])
  LOF[[t]] <- g
}
if (is.null(edges)) edges <- data.frame(GENE = character(0), PROGRAM = character(0), beta = numeric(0),
                                        padj = numeric(0), type = character(0), trait = character(0))
edges$key <- paste(edges$GENE, edges$PROGRAM, edges$type)
edges <- edges[!duplicated(edges$key), ]

## conservation counts
gcount <- table(unlist(lapply(conc, unique)))                    # gene -> #traits concordant
prognodes <- sort(unique(c(unlist(selP), unlist(selR))))
prognodes <- prognodes[order(as.numeric(sub("P", "", prognodes)))]
pcount <- sapply(prognodes, function(p) sum(sapply(TRAITS, function(t) p %in% c(selP[[t]], selR[[t]]))))
names(pcount) <- prognodes

## optional: keep only annotated program nodes. Filter here so mem2, the regulates edges, p2t and the
## gene set all inherit it; genes left with no edge to a shown program drop out in the layout step.
if (ANNOTATED_ONLY) {
  keepP <- prognodes[prognodes %in% names(labmap) & !is.na(labmap[prognodes]) & labmap[prognodes] != ""]
  prognodes <- keepP; pcount <- pcount[prognodes]
  edges <- edges[edges$PROGRAM %in% prognodes, , drop = FALSE]
}

## membership edges for every program node (union of concordant genes; gene-list rule unchanged)
conc_union <- unique(unlist(conc, use.names = FALSE))
mem2 <- do.call(rbind, lapply(prognodes, function(p) {
  tg <- names(sort(GEPs[, p], decreasing = TRUE))[1:TOP]
  hit <- intersect(conc_union, tg)
  if (length(hit)) data.frame(GENE = hit, PROGRAM = p, beta = NA_real_, padj = NA_real_, type = "member",
                              trait = NA_character_, stringsAsFactors = FALSE) else NULL
}))
if (!is.null(mem2) && nrow(mem2)) {
  mem2$key <- paste(mem2$GENE, mem2$PROGRAM, mem2$type)
  mem2 <- mem2[!(mem2$key %in% edges$key), , drop = FALSE]
  if (nrow(mem2)) edges <- rbind(edges, mem2[, colnames(edges), drop = FALSE])
}
mm <- edges$type == "member"
edges$strength <- abs(edges$beta)
edges$strength[mm] <- mapply(function(g, p) GEPs[g, p], edges$GENE[mm], edges$PROGRAM[mm])

## ---------------- program -> trait strength ----------------
p2t <- data.frame()
for (t in TRAITS) {
  regs <- selR[[t]]
  df <- merge(LOF[[t]], Breg, by = "GENE"); df <- merge(df, shet, by = "ensg")
  df <- df[, c("post_mean", regs, "shet")]
  for (i in 1:ncol(df)) { v <- df[, i]; v[is.infinite(v)] <- max(v[!is.infinite(v)]); df[, i] <- scale(v) }
  co <- summary(lm(post_mean ~ ., data = df))$coefficients
  for (p in regs) if (p %in% rownames(co))
    p2t <- rbind(p2t, data.frame(PROGRAM = p, trait = t, sign = sign(co[p, "Estimate"]),
                                 mlogp = -log10(max(co[p, "Pr(>|t|)"], 1e-300)), route = "regulator"))
  sh <- shet[is.element(shet$ensg, LOF[[t]]$ensg), ]
  for (p in setdiff(selP[[t]], regs)) {
    tmp <- data.frame(ensg = row.names(GEPe), GEPscore = GEPe[, p], stringsAsFactors = FALSE)
    d  <- merge(tmp, LOF[[t]], by = "ensg"); d2 <- d[order(d$GEPscore, decreasing = TRUE), ][1:TOP, ]
    obs <- mean(d2$post_mean)
    st <- sh[is.element(sh$ensg, d2$ensg), ]; sa <- sh[is.element(sh$ensg, d$ensg), ]
    A <- table(st$shet_BIN); rnd <- numeric(NPERM_SIGN)
    for (I in 1:NPERM_SIGN) { set.seed(I); gs <- c()
      for (q in seq_along(A)) gs <- c(gs, sample(sa$ensg[is.element(sa$shet_BIN, names(A)[q])], A[q]))
      rnd[I] <- mean(d[is.element(d$ensg, gs), "post_mean"]) }
    P1 <- (rank(c(obs, rnd))[1] / length(rnd)) * 2; P2 <- ((length(rnd) + 2 - rank(c(obs, rnd))[1]) / length(rnd)) * 2
    p2t <- rbind(p2t, data.frame(PROGRAM = p, trait = t, sign = sign(obs - mean(rnd)),
                                 mlogp = -log10(max(min(P1, P2, 1), 1 / NPERM_SIGN)), route = "program"))
  }
}
if (ANNOTATED_ONLY) p2t <- p2t[p2t$PROGRAM %in% prognodes, , drop = FALSE]

## ---------------- layout ----------------
genes <- sort(unique(edges$GENE[!is.na(edges$GENE)]))
if (length(genes) == 0 || length(prognodes) == 0) {
  cat(sprintf("[multi] %s | %s: nothing to draw after filtering (genes=%d, programs=%d) — skipped\n",
              CT, GROUP, length(genes), length(prognodes)))
  quit(save = "no", status = 0)
}
gc_of <- function(g) as.integer(gcount[g]); gc_of_v <- sapply(genes, function(g) if (g %in% names(gcount)) as.integer(gcount[g]) else 0L)
gclass <- ifelse(gc_of_v >= 3, "conserved (>=3 traits)", ifelse(gc_of_v == 2, "shared (2 traits)", "trait-specific"))
names(gclass) <- genes
ord <- order(-gc_of_v, genes); genes <- genes[ord]; gclass <- gclass[genes]
gy <- setNames(seq_along(genes), genes)
py <- setNames(seq(1, length(genes), length.out = length(prognodes)), prognodes)
## N trait nodes spread across the gene axis
nT <- length(TRAITS)
ty <- setNames(seq(length(genes) * 0.85, length(genes) * 0.15, length.out = nT), TRAITS)

rel <- function(v) { v <- abs(v); if (all(!is.finite(v)) || max(v, na.rm = TRUE) == 0) return(rep(0.5, length(v)))
                     pmax(v / max(v, na.rm = TRUE), 0.06) }
ed <- data.frame(x = 0, y = gy[edges$GENE], xend = 1.25, yend = py[edges$PROGRAM],
                 dir = ifelse(edges$type == "member", "membership",
                              ifelse(!is.na(edges$beta) & edges$beta > 0, "up-regulates", "down-regulates")),
                 stringsAsFactors = FALSE)
ed$rel <- NA_real_
ed$rel[edges$type == "member"]    <- rel(edges$strength[edges$type == "member"])    * 0.32
ed$rel[edges$type == "regulates"] <- rel(edges$strength[edges$type == "regulates"]) * 0.85
ed$alpha <- ifelse(edges$type == "member", 0.28, 0.60)
ed <- ed[!is.na(ed$y) & !is.na(ed$yend), ]

pt <- data.frame(x = 1.25, y = py[p2t$PROGRAM], xend = 2.7, yend = ty[p2t$trait],
                 dir = ifelse(p2t$sign > 0, "increases trait", "decreases trait"),
                 rel = 0.40 + 0.60 * rel(p2t$mlogp))

gn <- data.frame(x = 0, y = gy[genes], lab = genes, cls = gclass[genes],
                 big = gc_of_v[genes] >= 2)
pn <- data.frame(x = 1.25, y = py[prognodes],
                 lab = paste0(prognodes, ifelse(prognodes %in% names(labmap) & labmap[prognodes] != "",
                                                 paste0("  ", labmap[prognodes]), "")),
                 cls = ifelse(pcount[prognodes] >= 3, "conserved (>=3 traits)",
                       ifelse(pcount[prognodes] == 2, "shared (2 traits)", "one trait")))
tn <- data.frame(x = 2.7, y = ty, lab = sapply(TRAITS, short_of))

gene_fill <- c("conserved (>=3 traits)" = "#d1495b", "shared (2 traits)" = "#edae49", "trait-specific" = "grey90")
prog_fill <- c("conserved (>=3 traits)" = "#2e7d32", "shared (2 traits)" = "#7ac77a", "one trait" = "grey93")

g <- ggplot() + theme_void(base_size = 13) +
  geom_segment(data = ed, aes(x, y, xend = xend, yend = yend, color = dir, linewidth = rel, alpha = alpha)) +
  geom_segment(data = pt, aes(x, y, xend = xend, yend = yend, color = dir, linewidth = rel), alpha = 0.85,
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  scale_alpha_identity() +
  geom_point(data = gn, aes(x, y, fill = cls), shape = 21, size = 3.0, color = "grey25") +
  geom_text(data = gn, aes(x - 0.04, y, label = lab, fontface = ifelse(big, "bold", "plain"),
                           size = ifelse(big, 2.9, 2.2)), hjust = 1) +
  geom_label(data = pn, aes(x, y, label = lab, fill = cls), size = 2.9, label.size = 0.3) +
  geom_label(data = tn, aes(x, y, label = lab), size = 3.8, fill = "#ffe3e3", fontface = "bold") +
  scale_size_identity() +
  scale_linewidth_continuous(name = "relative strength (within edge class)", range = c(0.12, 2.4),
                             limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1.0)) +
  scale_color_manual(name = NULL,
    values = c("up-regulates" = "#c0392b", "down-regulates" = "#2c6fbb", "membership" = "grey72",
               "increases trait" = "#c0392b", "decreases trait" = "#2c6fbb")) +
  scale_fill_manual(name = NULL, values = c(gene_fill, prog_fill)) +
  guides(fill = guide_legend(override.aes = list(label = "", size = 4), order = 1, nrow = 2),
         color = guide_legend(override.aes = list(label = "", linewidth = 1.2), order = 2, nrow = 1),
         linewidth = guide_legend(order = 3, nrow = 1, override.aes = list(colour = "grey30", label = ""))) +
  coord_cartesian(xlim = c(-0.85, 3.2), clip = "off") +
  labs(title = sprintf("%s | %s | conserved programs & shared regulator genes across %d traits%s",
                       LABEL, GROUP, nT, if (ANNOTATED_ONLY) "  (annotated programs only)" else ""),
       subtitle = sprintf(paste0("traits: %s\n",
         "genes = union of concordant top hits (each trait's own |gamma|>q99 rule); colour = # traits concordant. ",
         "programs coloured by # traits selecting them.\n",
         "LINE WIDTH = strength within edge class: regulates |beta_x(P)|; program->trait -log10(P); ",
         "membership = cNMF loading (compressed)."),
         paste(sapply(TRAITS, short_of), collapse = ", "))) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.spacing.y = unit(1, "pt"),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(10, 12, 8, 34),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 8.3, colour = "grey30"))

h <- max(7.5, length(genes) * 0.23)
base <- sprintf("Fig5a_%s_%s%s", GROUP, CT, if (ANNOTATED_ONLY) "_annot" else "")
ggsave(file.path(FIGDIR, paste0(base, ".png")), g, width = 13.2, height = h, dpi = 150, limitsize = FALSE, bg = "white")
ggsave(file.path(FIGDIR, paste0(base, ".pdf")), g, width = 13.2, height = h, limitsize = FALSE, bg = "white")

## ---------------- companion tables: conserved programs + shared regulator genes ----------------
prog_tab <- data.frame(PROGRAM = prognodes, n_traits = as.integer(pcount[prognodes]),
                       label = ifelse(prognodes %in% names(labmap), labmap[prognodes], ""),
                       traits = sapply(prognodes, function(p)
                         paste(sapply(TRAITS[sapply(TRAITS, function(t) p %in% c(selP[[t]], selR[[t]]))], short_of), collapse = ";")),
                       stringsAsFactors = FALSE)
prog_tab <- prog_tab[order(-prog_tab$n_traits, prog_tab$PROGRAM), ]
write.table(prog_tab, file.path(FIGDIR, paste0(base, "_programs.tsv")), row.names = FALSE, sep = "\t", quote = FALSE)

gene_tab <- data.frame(GENE = names(gcount), n_traits = as.integer(gcount),
                       traits = sapply(names(gcount), function(gg)
                         paste(sapply(TRAITS[sapply(TRAITS, function(t) gg %in% conc[[t]])], short_of), collapse = ";")),
                       stringsAsFactors = FALSE)
gene_tab <- gene_tab[order(-gene_tab$n_traits, gene_tab$GENE), ]
write.table(gene_tab, file.path(FIGDIR, paste0(base, "_genes.tsv")), row.names = FALSE, sep = "\t", quote = FALSE)

cat(sprintf("[multi] %s | %s: %d traits, %d union genes (%d shared>=2), %d programs (%d shared>=2)\n",
            CT, GROUP, nT, length(genes), sum(gc_of_v >= 2), length(prognodes), sum(pcount >= 2)))
cat(sprintf("   conserved programs (>=2 traits): %s\n",
            paste(prog_tab$PROGRAM[prog_tab$n_traits >= 2], collapse = ", ")))
cat(sprintf("   shared regulator genes (>=2 traits): %s\n",
            paste(gene_tab$GENE[gene_tab$n_traits >= 2], collapse = ", ")))
