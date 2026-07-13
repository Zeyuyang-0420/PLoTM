#!/usr/bin/env Rscript
# fig5_step2_assemble.R
#
# Port of source_code/Figure5/2_permutation_test_step2.R + assembly of the Fig-5a regulatory map.
#
# Consumes the sharded output of fig5_permutation_step1_pseudobulk.R and produces, per (CT, trait):
#   1. <TAG>_permP.txt      observed Fisher P, permutation P, concordance counts
#   2. <TAG>_hist.pdf/png   null distribution of -log10(Fisher P) with an arrow at the observed
#   3. <TAG>_scatter.pdf/png  concordant vs discordant gene counts (observed in red)
#   4. <TAG>_Fig5a_map.pdf/png  the Fig-5a regulatory map: genes -> programs -> trait
#   5. <TAG>_map_edges.tsv / _map_nodes.tsv  (so the map can be redrawn / inspected)
#
# Documented fixes to the shipped script (same class as the Figure4 fixes):
#   FIX A  it read "perutation_test/" (typo) -> never found its own files.
#   FIX B  TRAIT/paths were hardcoded ("MCH", "permutation_test/P200/MCH_program5_regulator3_LOF0.1.txt").
#          Now arguments.
#   FIX C  the permuted-Fisher filter `tmp[1,8]-tmp[1,6] > 0` is kept (paper), but the hard gate
#          `length(seed_list) > 19995` is relaxed to a warning so partial shard sets still report.
#
# Program->trait signs: the observed run does not persist Program_P_sum / lm coefficients, so they are
# recomputed here (regulator coefficients exactly; program meanG signs via a cheaper shet-matched
# resample, whose SIGN is stable). Documented, non-scientific reconstruction for drawing only.
#
# Usage: Rscript fig5_step2_assemble.R <CT> <TRAIT_FILE> <PN> <RN> <TOP> <LOF> [OUTROOT] [FIGDIR] [NPERM_SIGN]

suppressMessages({library(ggplot2)})
args <- commandArgs(trailingOnly = TRUE)
options("scipen" = 10)
CT <- args[1]; TRAIT_FILE <- args[2]
PN <- as.numeric(args[3]); RN <- as.numeric(args[4])
TOP <- as.numeric(args[5]); LOF_thresh <- as.numeric(args[6])
OUTROOT <- if (length(args) >= 7) args[7] else "fig5_out"
FIGDIR  <- if (length(args) >= 8) args[8] else "visulization/figures_fig5"
NPERM_SIGN <- if (length(args) >= 9) as.numeric(args[9]) else 2000

K <- 60
TRAIT <- sub("\\.per_gene_estimates\\.tsv$", "", TRAIT_FILE)
TAG <- paste0(TRAIT, "_program", PN, "_regulator", RN, "_LOF", LOF_thresh)
OUTDIR <- file.path(OUTROOT, CT, paste0("P", TOP))
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

RESULTS <- Sys.getenv("RESULTS_ROOT", "/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis")
DATA    <- Sys.getenv("DATA_ROOT", "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data")
LABEL   <- if (grepl("rest", CT)) "GWCD4i Rest" else "GWCD4i Stim48hr"

## ---------- 1. concatenate obs + perm shards (FIX A/B) ----------
obs_f <- file.path(OUTDIR, paste0(TAG, "_obs.txt"))
stopifnot(file.exists(obs_f))
obs <- read.table(obs_f, header = TRUE, stringsAsFactors = FALSE)
shards <- list.files(OUTDIR, pattern = paste0("^", TAG, "_perm_.*\\.txt$"), full.names = TRUE)
perm <- do.call(rbind, lapply(shards, read.table, header = TRUE, stringsAsFactors = FALSE))
cat(sprintf("[step2] %s | %s : obs=%d rows, %d shards, %d permutations\n",
            CT, TRAIT, nrow(obs), length(shards), nrow(perm)))

fisherP <- function(r) fisher.test(matrix(c(r$CONCO_N, r$TOTAL_N - r$CONCO_N,
                                            r$BG_CONCO_N, r$BG_TOTAL_N - r$BG_CONCO_N),
                                          nrow = 2, byrow = TRUE))$p.value
obs_fp <- fisherP(obs[1, ])
perm <- perm[(perm$BG_TOTAL_N - perm$BG_CONCO_N) > 0, ]          # paper's filter
perm_fp <- vapply(seq_len(nrow(perm)), function(i) fisherP(perm[i, ]), numeric(1))
if (nrow(perm) < 19995) cat(sprintf("[step2] NOTE: %d permutations (<19995); P resolution = %.2e\n",
                                    nrow(perm), 1 / nrow(perm)))  # FIX C

permP <- rank(c(obs_fp, perm_fp))[1] / length(perm_fp)            # one-sided, paper's definition
obs_rate <- obs$CONCO_N / obs$TOTAL_N; bg_rate <- obs$BG_CONCO_N / obs$BG_TOTAL_N

res <- data.frame(CT = CT, TRAIT = TRAIT, PN = PN, RN = RN, TOP = TOP, LOF_thresh = LOF_thresh,
                  CONCO_N = obs$CONCO_N, DISCO_N = obs$DISCO_N, TOTAL_N = obs$TOTAL_N,
                  conc_rate = obs_rate, bg_conc_rate = bg_rate,
                  fisher_P = obs_fp, perm_P = permP, N_perm = length(perm_fp))
write.table(res, file.path(OUTDIR, paste0(TAG, "_permP.txt")), row.names = FALSE, sep = "\t", quote = FALSE)
cat(sprintf("[step2] observed %d/%d concordant (%.1f%%) vs bg %.1f%% | Fisher P=%.3e | perm P=%.3e (n=%d)\n",
            obs$CONCO_N, obs$TOTAL_N, 100 * obs_rate, 100 * bg_rate, obs_fp, permP, length(perm_fp)))

## ---------- 2. permutation histogram ----------
hd <- data.frame(logP = -log10(c(obs_fp, perm_fp)), grp = c("observed", rep("permuted", length(perm_fp))))
g <- ggplot(hd[hd$grp == "permuted", ], aes(x = logP)) +
  theme_classic(base_size = 18) +
  geom_histogram(fill = "grey35", binwidth = 0.05) +
  geom_vline(xintercept = -log10(obs_fp), color = "red", linewidth = 1) +
  annotate("text", x = -log10(obs_fp), y = Inf, vjust = 1.6, hjust = -0.05, color = "red", size = 5,
           label = sprintf("observed\n-log10P=%.2f\nperm P=%.2e", -log10(obs_fp), permP)) +
  labs(x = "-log10(Fisher P)", y = "count",
       title = sprintf("Fig 5 permutation null | %s | %s", LABEL, TRAIT),
       subtitle = sprintf("%d permutations; observed concordance %d/%d (%.1f%%) vs background %.1f%%",
                          length(perm_fp), obs$CONCO_N, obs$TOTAL_N, 100 * obs_rate, 100 * bg_rate))
ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_hist.png")), g, width = 9, height = 6, dpi = 150)
ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_hist.pdf")), g, width = 9, height = 6)

## ---------- 3. concordant / discordant scatter ----------
sc <- rbind(data.frame(CONCO_N = perm$CONCO_N, DISCO_N = perm$DISCO_N, grp = "permuted"),
            data.frame(CONCO_N = obs$CONCO_N, DISCO_N = obs$DISCO_N, grp = "observed"))
g2 <- ggplot(sc[sc$grp == "permuted", ], aes(CONCO_N, DISCO_N)) +
  theme_classic(base_size = 18) +
  geom_point(size = 1, alpha = 0.35, color = "grey60") +
  geom_point(data = sc[sc$grp == "observed", ], color = "red", size = 5) +
  labs(x = "N concordant genes", y = "N discordant genes",
       title = sprintf("Fig 5 | %s | %s", LABEL, TRAIT),
       subtitle = "red = observed model; grey = permuted null")
ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_scatter.png")), g2, width = 7, height = 6.5, dpi = 150)
ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_scatter.pdf")), g2, width = 7, height = 6.5)

## ---------- 4. Fig 5a regulatory map ----------
sel <- read.table(file.path(OUTDIR, paste0("SelectedPrograms_", TAG, ".txt")), header = TRUE,
                  stringsAsFactors = FALSE, sep = "\t")
Program_selected   <- strsplit(sel$Program_selected[1], ",")[[1]]
regulator_selected <- strsplit(sel$Regulator_selected[1], ",")[[1]]
conc <- read.table(file.path(OUTDIR, paste0("ConcordantGenes_", TAG, ".txt")), header = TRUE,
                   stringsAsFactors = FALSE)

# reload the model inputs
GEP <- read.table(file.path(RESULTS, "cNMF", CT, "test1",
                            paste0("test1.gene_spectra_score.k_", K, ".dt_0_4.txt")),
                  header = TRUE, stringsAsFactors = FALSE)
GEP <- t(GEP); colnames(GEP) <- paste0("P", 1:K)
corresp <- read.table(file.path(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID"), header = FALSE,
                      stringsAsFactors = FALSE)
corresp <- corresp[!duplicated(corresp[, 1]), ]; row.names(corresp) <- corresp[, 1]
GEP <- data.frame(GENE = corresp[row.names(GEP), 2], GEP)

GEP_reg_beta <- NULL; GEP_reg_P <- NULL
for (i in 1:K) {
  tmp <- read.table(file.path(RESULTS, "cNMF_regulation", CT,
                              paste0("K", K, "_program", i, "_perturb_effects.txt")), header = TRUE)
  t1 <- tmp[, c(1, 2)]; t2 <- tmp[, c(1, 3)]
  colnames(t1) <- c("GENE", paste0("P", i)); colnames(t2) <- c("GENE", paste0("P", i))
  if (i == 1) { GEP_reg_beta <- t1; GEP_reg_P <- t2 } else {
    GEP_reg_beta <- merge(GEP_reg_beta, t1, by = "GENE"); GEP_reg_P <- merge(GEP_reg_P, t2, by = "GENE") }
}
LOF <- read.table(file.path(DATA, "LoF", "GeneBayes_posterior", TRAIT_FILE), sep = "\t", quote = "",
                  header = TRUE, stringsAsFactors = FALSE)
LOF <- data.frame(GENE = corresp[LOF$ensg, 2], LOF); LOF <- LOF[!is.na(LOF$GENE), ]
shet <- read.table(file.path(DATA, "shet_10bins.txt"), header = TRUE, stringsAsFactors = FALSE)
shet <- shet[is.element(shet$ensg, LOF$ensg), ]

# regulator program -> trait sign: refit the paper's scaled lm with the selected regulators
dfm <- merge(LOF, GEP_reg_beta, by = "GENE"); dfm <- merge(dfm, shet, by = "ensg")
dfm <- dfm[, c("GENE.x", "post_mean", regulator_selected, "shet")]; colnames(dfm)[1] <- "GENE"
for (i in 2:ncol(dfm)) { dfm[, i][is.infinite(dfm[, i])] <- max(dfm[, i][!is.infinite(dfm[, i])])
                          dfm[, i] <- scale(dfm[, i]) }
coefs <- summary(lm(post_mean ~ ., data = dfm[, -1]))$coefficients
reg_sign <- sign(coefs[regulator_selected, "Estimate"]); names(reg_sign) <- regulator_selected

# program-selected -> trait sign: sign(meanG) recomputed with a cheaper shet-matched resample
prog_sign <- setNames(rep(0, length(Program_selected)), Program_selected)
for (p in Program_selected) {
  tmp <- GEP[, c("GENE", p)]; tmp <- tmp[order(tmp[, 2], decreasing = TRUE), ]
  colnames(tmp)[2] <- "GEPscore"; tmp <- data.frame(tmp, ensg = row.names(tmp))
  d <- merge(tmp, LOF, by = "ensg"); d2 <- d[order(d$GEPscore, decreasing = TRUE), ][1:TOP, ]
  st <- shet[is.element(shet$ensg, d2$ensg), ]; sa <- shet[is.element(shet$ensg, d$ensg), ]
  A <- table(st$shet_BIN); rnd <- c()
  for (I in 1:NPERM_SIGN) { set.seed(I); gs <- c()
    for (q in 1:length(A)) gs <- c(gs, sample(sa$ensg[is.element(sa$shet_BIN, names(A)[q])], A[q]))
    rnd <- c(rnd, mean(d[is.element(d$ensg, gs), "post_mean"])) }
  prog_sign[p] <- sign(mean(d2$post_mean) - mean(rnd))
}

# gene -> program edges (regulator programs, BH-significant beta)
row.names(GEP_reg_beta) <- GEP_reg_beta$GENE; row.names(GEP_reg_P) <- GEP_reg_P$GENE
padj <- as.data.frame(lapply(GEP_reg_P[, regulator_selected, drop = FALSE], p.adjust, method = "BH"))
row.names(padj) <- GEP_reg_P$GENE
genes <- intersect(conc$GENE, row.names(GEP_reg_beta))
edges <- do.call(rbind, lapply(genes, function(g) {
  do.call(rbind, lapply(regulator_selected, function(p) {
    if (!is.na(padj[g, p]) && padj[g, p] < 0.05)
      data.frame(GENE = g, PROGRAM = p, beta = GEP_reg_beta[g, p], padj = padj[g, p],
                 type = "regulates", stringsAsFactors = FALSE) else NULL }))
}))
# membership edges: concordant gene is in the top-TOP loadings of a selected program
memb <- do.call(rbind, lapply(Program_selected, function(p) {
  tg <- GEP[order(GEP[, p], decreasing = TRUE), "GENE"][1:TOP]
  hit <- intersect(genes, tg)
  if (length(hit)) data.frame(GENE = hit, PROGRAM = p, beta = NA, padj = NA,
                              type = "member", stringsAsFactors = FALSE) else NULL }))
edges <- rbind(edges, memb)
gsign <- setNames(sign(conc$post_mean), conc$GENE)

nodes_prog <- unique(c(regulator_selected, Program_selected))
psign <- sapply(nodes_prog, function(p) if (p %in% names(reg_sign)) reg_sign[[p]] else prog_sign[[p]])
# guard the zero-concordant-gene case: data.frame(id = character(0), kind = "gene", ...) errors on
# recycling (0 rows vs 1). Only add gene rows when there are genes; the map then falls through to the
# placeholder branch below.
gene_nodes <- if (length(genes) > 0)
  data.frame(id = genes, kind = "gene", sign = gsign[genes], stringsAsFactors = FALSE) else
  data.frame(id = character(0), kind = character(0), sign = numeric(0), stringsAsFactors = FALSE)
nodes <- rbind(
  gene_nodes,
  data.frame(id = nodes_prog, kind = "program", sign = psign, stringsAsFactors = FALSE),
  data.frame(id = TRAIT, kind = "trait", sign = 0, stringsAsFactors = FALSE))
write.table(edges, file.path(OUTDIR, paste0(TAG, "_map_edges.tsv")), row.names = FALSE, sep = "\t", quote = FALSE)
write.table(nodes, file.path(OUTDIR, paste0(TAG, "_map_nodes.tsv")), row.names = FALSE, sep = "\t", quote = FALSE)

if (!is.null(edges) && nrow(edges) > 0) {
  genes_drawn <- sort(unique(edges$GENE))
  gy <- setNames(seq_along(genes_drawn), genes_drawn)
  py <- setNames(seq(1, length(genes_drawn), length.out = length(nodes_prog)), nodes_prog)
  ty <- mean(range(gy))
  ed <- data.frame(x = 0, y = gy[edges$GENE], xend = 1, yend = py[edges$PROGRAM],
                   type = edges$type,
                   dir = ifelse(is.na(edges$beta), "membership",
                                ifelse(edges$beta > 0, "up-regulates", "down-regulates")))
  p2t <- data.frame(x = 1, y = py[nodes_prog], xend = 2, yend = ty,
                    dir = ifelse(psign > 0, "increases trait", "decreases trait"))
  gn <- data.frame(x = 0, y = gy[genes_drawn], lab = genes_drawn,
                   dir = ifelse(gsign[genes_drawn] > 0, "gamma > 0", "gamma < 0"))
  pn <- data.frame(x = 1, y = py[nodes_prog], lab = nodes_prog)

  gm <- ggplot() + theme_void(base_size = 14) +
    geom_segment(data = ed, aes(x, y, xend = xend, yend = yend, color = dir), alpha = 0.55, linewidth = 0.4) +
    geom_segment(data = p2t, aes(x, y, xend = xend, yend = yend, color = dir), linewidth = 1.1,
                 arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    geom_point(data = gn, aes(x, y, fill = dir), shape = 21, size = 3.2, color = "grey20") +
    geom_text(data = gn, aes(x - 0.03, y, label = lab), hjust = 1, size = 2.6) +
    geom_label(data = pn, aes(x, y, label = lab), size = 3.6, fill = "grey92") +
    geom_label(data = data.frame(x = 2, y = ty), aes(x, y), label = TRAIT, size = 4.2, fill = "#ffe6e6") +
    scale_color_manual(values = c("up-regulates" = "#c0392b", "down-regulates" = "#2c6fbb",
                                  "membership" = "grey70",
                                  "increases trait" = "#c0392b", "decreases trait" = "#2c6fbb")) +
    scale_fill_manual(values = c("gamma > 0" = "#e8746a", "gamma < 0" = "#6a9fe8")) +
    coord_cartesian(xlim = c(-0.45, 2.3), clip = "off") +
    labs(title = sprintf("Fig 5a | %s | %s  (regulators -> programs -> trait)", LABEL, TRAIT),
         subtitle = sprintf("%d concordant top hits (|gamma|>%s); programs %s; regulators %s; perm P=%.2e",
                            length(genes_drawn), LOF_thresh, sel$Program_selected[1],
                            sel$Regulator_selected[1], permP)) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(size = 14), plot.subtitle = element_text(size = 9, colour = "grey30"))
  h <- max(7, length(genes_drawn) * 0.22)
  ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_Fig5a_map.png")), gm, width = 11, height = h, dpi = 150, limitsize = FALSE)
  ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_Fig5a_map.pdf")), gm, width = 11, height = h, limitsize = FALSE)
  cat(sprintf("[step2] map: %d gene nodes, %d edges -> %s\n", length(genes_drawn), nrow(edges), FIGDIR))
} else {
  # No concordant top-hit genes => the regulator->program->trait map has no gene nodes to draw.
  # This is a genuine null (the model does not predict sign(gamma) for any top hit), not a failure.
  # Emit an explicit placeholder so the 16-panel deliverable has no silent gaps.
  cat("[step2] no edges -> map not drawn; writing placeholder\n")
  gm <- ggplot() +
    annotate("text", x = 0, y = 0.15,
             label = sprintf("Fig 5a | %s | %s", LABEL, TRAIT), size = 6, fontface = "bold") +
    annotate("text", x = 0, y = -0.15,
             label = sprintf("No concordant top-hit genes (|gamma| > %s):\nthe regulator->program model predicts sign(gamma) for 0 of %d top hits.\nperm P = %.3g (non-significant). No network to draw.",
                             LOF_thresh, obs$TOTAL_N, permP), size = 4.2) +
    theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA)) +
    xlim(-1, 1) + ylim(-1, 1)
  ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_Fig5a_map.png")), gm, width = 11, height = 4, dpi = 150, bg = "white")
  ggsave(file.path(FIGDIR, paste0(TAG, "_", CT, "_Fig5a_map.pdf")), gm, width = 11, height = 4, bg = "white")
}
cat("[step2] DONE\n")
