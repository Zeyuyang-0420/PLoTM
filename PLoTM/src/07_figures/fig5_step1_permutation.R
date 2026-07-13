#!/usr/bin/env Rscript
# fig5_permutation_step1_pseudobulk.R
#
# Port of source_code/Figure5/1_permutation_test_step1.R to the GWCD4i pseudobulk cNMF runs and the
# RA LoF traits. Builds the Fig-5a "regulators -> programs -> trait" model:
#   step1  select top Program_N programs by program-burden (shet-matched permutation)
#   step2  select top Regulator_N programs by regulator-burden via leaps::regsubsets (best subset)
#   step3  predict sign(gamma) per gene (program membership overrides regulator effect) and count
#          concordant / discordant genes among top hits (|gamma| > LOF_thresh) vs background
# then repeats the whole thing on NPERM label-permuted gamma vectors to build the null.
#
# ALL STATISTICS ARE THE PAPER'S. Only the following are changed, each documented:
#   FIX 1  `dir.create(paste0(..., Program_top_def))` was called on line 5 BEFORE Program_top_def was
#          assigned (line 9) -> "object not found". Moved after the args are parsed.
#   FIX 2  LOF was built as data.frame(gene=...) (lowercase) but then merged `by="GENE"` -> error.
#          The column is now named GENE.
#   FIX 3  the IRF branch assigned an undefined object `Ret`; TRAIT was hardcoded "MCH". TRAIT and the
#          gamma file are now arguments (no dead branches).
#   FIX 4  input paths hardcoded K562GW / cNMF_all...dt_0_5. Now parameterised: our essential-style
#          consensus is dt_0_4 and CT-named.
#   FIX 5  the permutation loop re-wrote the whole summary table to disk on EVERY iteration. Now
#          written every 250 iterations and at the end (pure I/O; no statistic changed).
#
# HARNESS PARALLELIZATION (non-scientific): MODE splits the run into the expensive OBSERVED pass and
# independent PERMUTATION shards. Every permutation calls set.seed(SEED) and re-selects its own
# programs/regulators from the shuffled gamma, so sharding seed ranges across processes is exactly
# equivalent to the serial 1..NPERM loop. Shards are concatenated afterwards for step 2.
#
# Usage:
#   Rscript fig5_permutation_step1_pseudobulk.R <CT> <TRAIT_FILE> <Program_N> <Regulator_N> \
#           <Program_top_def> <LOF_thresh> <MODE:obs|perm> [SEED_START] [SEED_END] \
#           [NPERM_OBS=10000] [OUTROOT=permutation_test]

args <- commandArgs(trailingOnly = TRUE)
options("scipen" = 10)

CT              <- as.character(args[1])              # e.g. GWCD4i_rest_pseudobulk
TRAIT_FILE      <- as.character(args[2])              # e.g. Genebass_RA_M06.per_gene_estimates.tsv
Program_N       <- as.numeric(args[3])                # 1-6
Regulator_N     <- as.numeric(args[4])                # 1-6
Program_top_def <- as.numeric(args[5])                # 100, 200, 300
LOF_thresh      <- as.numeric(args[6])                # RA: 0.03 (paper used 0.1 for MCH; see README)
MODE            <- if (length(args) >= 7) as.character(args[7]) else "obs"   # "obs" | "perm"
SEED_START      <- if (length(args) >= 8) as.numeric(args[8]) else 1
SEED_END        <- if (length(args) >= 9) as.numeric(args[9]) else 0
NPERM_OBS       <- if (length(args) >= 10) as.numeric(args[10]) else 10000
OUTROOT         <- if (length(args) >= 11) as.character(args[11]) else "permutation_test"
stopifnot(MODE %in% c("obs", "perm"))

K     <- 60
TRAIT <- sub("\\.per_gene_estimates\\.tsv$", "", TRAIT_FILE)

RESULTS <- Sys.getenv("RESULTS_ROOT",
             "/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis")
DATA    <- Sys.getenv("DATA_ROOT",
             "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/data")

## FIX 1: dir.create AFTER args are parsed
OUTDIR <- file.path(OUTROOT, CT, paste0("P", Program_top_def))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("[fig5] CT=%s TRAIT=%s PN=%d RN=%d top=%d LOF=%s MODE=%s seeds=%d..%d NPERM_OBS=%d\n",
            CT, TRAIT, Program_N, Regulator_N, Program_top_def, LOF_thresh, MODE,
            SEED_START, SEED_END, NPERM_OBS))

## ---------------- data loading (FIX 4: CT-named, dt_0_4) ----------------
GEP <- read.table(file.path(RESULTS, "cNMF", CT, "test1",
                            paste0("test1.gene_spectra_score.k_", K, ".dt_0_4.txt")),
                  header = TRUE, stringsAsFactors = FALSE)
GEP <- t(GEP)
colnames(GEP) <- paste0("P", 1:K)
corresp <- read.table(file.path(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID"),
                      header = FALSE, stringsAsFactors = FALSE)
corresp <- corresp[!duplicated(corresp[, 1]), ]
row.names(corresp) <- corresp[, 1]
GEP <- data.frame(GENE = corresp[row.names(GEP), 2], GEP)

## ESTIMABILITY (no-op for pure Perturb-seq): the regulator candidate programs are the ESTIMABLE
## ones. With the disease-single-cell basis some programs have empty beta files; restrict the beta
## merge + regsubsets candidates to the file named by env ALIVE_PROGRAMS. Absent -> all 1:K (pure
## Perturb-seq unchanged). Program SELECTION (program-burden) still ranges over all K.
.alivef <- Sys.getenv("ALIVE_PROGRAMS")
if (!nzchar(.alivef)) for (cand in c(file.path(RESULTS, "cNMF", CT, "alive_programs.txt"),
                                     file.path("data/Perturbseq/cNMF", CT, "alive_programs.txt")))
  if (file.exists(cand)) { .alivef <- cand; break }
PROGS <- if (nzchar(.alivef) && file.exists(.alivef)) as.numeric(readLines(.alivef)) else 1:K
CAND_P <- paste0("P", PROGS)

GEP_reg_beta <- data.frame(); GEP_reg_P <- data.frame()
.first <- TRUE
for (i in PROGS) {
  tmp <- read.table(file.path(RESULTS, "cNMF_regulation", CT,
                              paste0("K", K, "_program", i, "_perturb_effects.txt")), header = TRUE)
  tmp1 <- tmp[, c(1, 2)]; tmp2 <- tmp[, c(1, 3)]
  colnames(tmp1) <- c("GENE", paste0("P", i))
  colnames(tmp2) <- c("GENE", paste0("P", i))
  if (.first) { GEP_reg_beta <- tmp1; GEP_reg_P <- tmp2; .first <- FALSE
  } else { GEP_reg_beta <- merge(GEP_reg_beta, tmp1, by = "GENE")
           GEP_reg_P    <- merge(GEP_reg_P,    tmp2, by = "GENE") }
}

## FIX 2 + FIX 3: gamma table keyed by GENE (upper case), trait chosen by argument
LOF_orig <- read.table(file.path(DATA, "LoF", "GeneBayes_posterior", TRAIT_FILE),
                       sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE)
LOF_orig <- data.frame(GENE = corresp[LOF_orig$ensg, 2], LOF_orig)
LOF_orig <- LOF_orig[!is.na(LOF_orig$GENE), ]

shet <- read.table(file.path(DATA, "shet_10bins.txt"), header = TRUE, stringsAsFactors = FALSE)
shet <- shet[is.element(shet$ensg, LOF_orig$ensg), ]

library(leaps)

## ---------------- shared helpers (verbatim statistics) ----------------
select_regulators <- function(LOF) {
  df <- merge(LOF, GEP_reg_beta, by = "GENE")
  df <- merge(df, shet, by = "ensg")
  df <- df[, c("post_mean", CAND_P, "shet")]
  for (i in 2:ncol(df)) df[, i][is.infinite(df[, i])] <- max(df[, i][!is.infinite(df[, i])])
  b <- leaps::regsubsets(post_mean ~ ., data = df, nbest = 1, nvmax = Regulator_N + 1, really.big = TRUE)
  w <- summary(b)[[1]]
  rs <- colnames(w)[w[nrow(w), ]]
  rs <- rs[!is.element(rs, c("(Intercept)", "shet"))]
  if (length(rs) != Regulator_N) {
    b <- leaps::regsubsets(post_mean ~ ., data = df, nbest = 1, nvmax = Regulator_N, really.big = TRUE)
    w <- summary(b)[[1]]
    rs <- colnames(w)[w[nrow(w), ]]
    rs <- rs[!is.element(rs, c("(Intercept)", "shet"))]
  }
  rs
}

score_model <- function(LOF, Program_selected, Program_P_sum, regulator_selected, SEEDLAB) {
  df <- merge(LOF, GEP_reg_beta, by = "GENE")
  df <- merge(df, shet, by = "ensg")
  df <- df[, c("GENE.x", "post_mean", regulator_selected, "shet")]
  colnames(df)[1] <- "GENE"
  for (i in 2:ncol(df)) {
    df[, i][is.infinite(df[, i])] <- max(df[, i][!is.infinite(df[, i])])
    df[, i] <- scale(df[, i])
  }
  fit1 <- lm(post_mean ~ ., data = df[, -1]); res1 <- summary(fit1)$coefficients

  df2 <- merge(df[, "GENE", drop = FALSE], GEP_reg_P, by = "GENE")
  row.names(df2) <- df2$GENE; df2 <- df2[df$GENE, ]

  sum_effect <- data.frame(GENE = df$GENE, effect = 0, regeffect = 0, proeffect = 0)
  for (p in regulator_selected) {
    tmp <- res1[p, "Estimate"] * df[, p]
    tmp[p.adjust(df2[, p], method = "BH") > 0.05] <- 0
    sum_effect$effect <- sum_effect$effect + tmp
  }
  sum_effect$regeffect <- sign(sum_effect$effect)
  for (p in rev(Program_selected)) {
    tmp  <- GEP[order(GEP[, p], decreasing = TRUE), "GENE"][1:Program_top_def]
    tmp2 <- sign(Program_P_sum$meanG[Program_P_sum$Program == gsub("P", "", p)])
    sum_effect$proeffect[is.element(sum_effect$GENE, tmp)] <- tmp2
  }

  LOF2 <- data.frame(LOF, BIN = ifelse(abs(LOF$post_mean) > LOF_thresh, 1, 2))
  BIN <- 1
  topgenes <- LOF2$GENE[LOF2$BIN == BIN]

  dff <- merge(LOF2[is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")], sum_effect, by = "GENE")
  dff <- data.frame(dff, post_mean_sign = sign(dff$post_mean))
  dff <- data.frame(dff, predicted_sign = dff$proeffect)
  dff$predicted_sign[dff$predicted_sign == 0] <- dff$regeffect[dff$predicted_sign == 0]
  CONCO_N <- length(which(dff$predicted_sign == dff$post_mean_sign))
  DISCO_N <- length(which(dff$predicted_sign == -dff$post_mean_sign & abs(dff$predicted_sign) > 0))
  TOTAL_N <- nrow(dff)

  dff2 <- merge(LOF2[!is.element(LOF2$GENE, topgenes), c("GENE", "post_mean")], sum_effect, by = "GENE")
  dff2 <- data.frame(dff2, post_mean_sign = sign(dff2$post_mean))
  dff2 <- data.frame(dff2, predicted_sign = dff2$proeffect)
  dff2$predicted_sign[dff2$predicted_sign == 0] <- dff2$regeffect[dff2$predicted_sign == 0]
  BG_CONCO_N <- length(which(dff2$predicted_sign == dff2$post_mean_sign))
  BG_DISCO_N <- length(which(dff2$predicted_sign == -dff2$post_mean_sign & abs(dff2$predicted_sign) > 0))
  BG_TOTAL_N <- nrow(dff2)

  list(row = data.frame(SEED = SEEDLAB, BIN = BIN, CONCO_N = CONCO_N, DISCO_N = DISCO_N,
                        TOTAL_N = TOTAL_N, BG_CONCO_N = BG_CONCO_N, BG_DISCO_N = BG_DISCO_N,
                        BG_TOTAL_N = BG_TOTAL_N),
       concordant = dff[dff$predicted_sign == dff$post_mean_sign, ])
}

summary_tab <- data.frame()
TAG <- paste0(TRAIT, "_program", Program_N, "_regulator", Regulator_N, "_LOF", LOF_thresh)

## ================= OBSERVED =================
if (MODE == "obs") {
LOF <- LOF_orig
Program_P_sum <- data.frame()
for (Program in 1:K) {
  tmp <- GEP[, c("GENE", paste0("P", Program))]
  tmp <- tmp[order(tmp[, 2], decreasing = TRUE), ]
  colnames(tmp)[2] <- "GEPscore"
  tmp <- data.frame(tmp, ensg = row.names(tmp))
  df  <- merge(tmp, LOF, by = "ensg")
  df2 <- df[order(df$GEPscore, decreasing = TRUE), ][1:Program_top_def, ]
  MEANgamma_top100 <- mean(df2[, "post_mean"])

  shet_tmp     <- shet[is.element(shet$ensg, df2$ensg), ]
  shet_tmp_all <- shet[is.element(shet$ensg, df$ensg), ]
  A <- table(shet_tmp$shet_BIN)
  random <- c()
  for (I in 1:NPERM_OBS) {
    set.seed(I)
    genes <- c()
    for (p in 1:length(A))
      genes <- c(genes, sample(shet_tmp_all$ensg[is.element(shet_tmp_all$shet_BIN, names(A)[p])], A[p]))
    random <- c(random, mean(df[is.element(df$ensg, genes), "post_mean"], na.omit = TRUE))
  }
  random_mean <- mean(random)
  P1 <- (rank(c(MEANgamma_top100, random))[1] / length(random)) * 2
  P2 <- ((length(random) + 2 - rank(c(MEANgamma_top100, random))[1]) / length(random)) * 2
  Program_P <- min(P1, P2)
  Program_P_sum <- rbind(Program_P_sum,
                         data.frame(Program = Program, P = Program_P, meanG = MEANgamma_top100 - random_mean))
}
Program_P_sum <- Program_P_sum[order(abs(Program_P_sum$meanG), decreasing = TRUE), ]
Program_P_sum <- Program_P_sum[order(Program_P_sum$P), ]
Program_selected <- paste0("P", Program_P_sum$Program[1:Program_N])

regulator_selected <- select_regulators(LOF)
obs <- score_model(LOF, Program_selected, Program_P_sum, regulator_selected, "TRUE")
summary_tab <- rbind(summary_tab, obs$row)
cat(sprintf("[fig5] OBSERVED programs=%s regulators=%s | CONCO=%d DISCO=%d TOTAL=%d (bg %d/%d/%d)\n",
            paste(Program_selected, collapse = ","), paste(regulator_selected, collapse = ","),
            obs$row$CONCO_N, obs$row$DISCO_N, obs$row$TOTAL_N,
            obs$row$BG_CONCO_N, obs$row$BG_DISCO_N, obs$row$BG_TOTAL_N))

write.table(data.frame(Program_selected = paste(Program_selected, collapse = ","),
                       Regulator_selected = paste(regulator_selected, collapse = ",")),
            file.path(OUTDIR, paste0("SelectedPrograms_", TAG, ".txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
# concordant top-hit genes -> the nodes drawn on the Fig 5a map
write.table(obs$concordant, file.path(OUTDIR, paste0("ConcordantGenes_", TAG, ".txt")),
            row.names = FALSE, sep = "\t", quote = FALSE)
OUTFILE <- file.path(OUTDIR, paste0(TAG, "_obs.txt"))
write.table(summary_tab, OUTFILE, row.names = FALSE, sep = "\t", quote = FALSE)
cat(sprintf("[fig5] DONE (obs) -> %s\n", OUTFILE))
}

## ================= PERMUTATIONS (sharded) =================
if (MODE == "perm") {
stopifnot(SEED_END >= SEED_START)
OUTFILE <- file.path(OUTDIR, sprintf("%s_perm_%d_%d.txt", TAG, SEED_START, SEED_END))
for (SEED in SEED_START:SEED_END) {
  LOF <- LOF_orig
  set.seed(SEED)
  LOF$post_mean <- sample(LOF$post_mean, nrow(LOF), replace = FALSE)

  Program_P_sum <- data.frame()
  for (Program in 1:K) {
    tmp <- GEP[, c("GENE", paste0("P", Program))]
    tmp <- tmp[order(tmp[, 2], decreasing = TRUE), ]
    colnames(tmp)[2] <- "GEPscore"
    tmp <- data.frame(tmp, ensg = row.names(tmp))
    df  <- merge(tmp, LOF, by = "ensg")
    o   <- order(abs(df$GEPscore), decreasing = TRUE)
    MEANgamma_top100 <- mean(df[o, "post_mean"][1:Program_top_def])
    Program_P <- wilcox.test(df[o, "post_mean"][1:Program_top_def],
                             df[o, "post_mean"][(Program_top_def + 1):nrow(df)])$p.value
    Program_P_sum <- rbind(Program_P_sum,
                           data.frame(Program = Program, P = Program_P,
                                      meanG = MEANgamma_top100 - mean(df$post_mean)))
  }
  Program_P_sum <- Program_P_sum[order(abs(Program_P_sum$meanG), decreasing = TRUE), ]
  Program_P_sum <- Program_P_sum[order(Program_P_sum$P), ]
  Program_selected <- paste0("P", Program_P_sum$Program[1:Program_N])

  regulator_selected <- select_regulators(LOF)
  summary_tab <- rbind(summary_tab,
                       score_model(LOF, Program_selected, Program_P_sum, regulator_selected, SEED)$row)

  ## FIX 5: write periodically instead of every iteration
  if (SEED %% 250 == 0) {
    write.table(summary_tab, OUTFILE, row.names = FALSE, sep = "\t", quote = FALSE)
    cat(sprintf("[fig5] perm %d (shard %d..%d)\n", SEED, SEED_START, SEED_END))
  }
}
write.table(summary_tab, OUTFILE, row.names = FALSE, sep = "\t", quote = FALSE)
cat(sprintf("[fig5] DONE (perm shard %d..%d) -> %s (%d rows)\n",
            SEED_START, SEED_END, OUTFILE, nrow(summary_tab)))
}
