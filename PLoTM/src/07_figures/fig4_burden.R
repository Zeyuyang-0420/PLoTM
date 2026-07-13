#!/usr/bin/env Rscript
# plot_fig4_burden.R
# Reproduce the paper's Figure 4C (Program burden effect vs Regulator-burden correlation) for our
# GWCD4i Stim48hr cNMF programs against RA rheumatoid-factor (Genebass_RF) loss-of-function gamma.
#
# Adapted from source_code/Figure4/2_plot_burden_program_regulators.R. The two axis scores are
# computed EXACTLY as the paper (signed -log10 P for program-burden and regulator-burden), and the
# paper's colour scheme is kept. Because RA rare-LoF burden is underpowered, NO program reaches
# Bonferroni significance (0.05/60); per request we therefore colour & label programs at nominal
# p < 0.05 and annotate the Bonferroni result explicitly. Both threshold lines are drawn (red dotted
# = Bonferroni, grey dashed = nominal p<0.05) so the significance context is unambiguous.

args <- commandArgs(trailingOnly = TRUE)
options("scipen" = 10)
FILE <- if (length(args) >= 1) args[1] else "Genebass_RF.per_gene_estimates.tsv"
K    <- if (length(args) >= 2) as.numeric(args[2]) else 60
DIR  <- if (length(args) >= 3) args[3] else
  "/mnt/scratch/ZY2/Hackathon/tcell_perturbseq/cNMF_RA_analysis/trait_association/GWCD4i_stim48_pseudobulk/ProgramLevel/Genebass_RF"
OUT  <- if (length(args) >= 4) args[4] else
  "/mnt/scratch/ZY2/Hackathon/Hackthron-Claude/Paper/reproduction/visulization"
LABEL <- if (length(args) >= 5) args[5] else "GWCD4i Stim48hr"   # condition label for title/filename
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

trait      <- sub("\\.per_gene_estimates\\.tsv$", "", FILE)       # e.g. Genebass_RF
label_slug <- gsub("[^A-Za-z0-9]+", "_", LABEL)                    # e.g. GWCD4i_Stim48hr

pro <- read.table(paste0(DIR, "/programs_enrichment_K",  K, "_", FILE), header = TRUE, stringsAsFactors = FALSE)
reg <- read.table(paste0(DIR, "/regulators_enrichment_K", K, "_", FILE), header = TRUE, stringsAsFactors = FALSE)
df  <- merge(pro, reg, by = "Program")

## ---- paper's signed -log10(P) scores (identical to Figure4/2_plot_burden_program_regulators.R) ----
df$Program_score   <- sign(df$MEANgamma_top100 - df$shet_adjusted_random_mean) * (-log10(df$MEANgamma_top100_shet_adjusted_P))
df$regulator_score <- sign(df$beta_withShet) * (-log10(df$P_withShet))

N    <- nrow(df)
bonf <- 0.05 / N        # Bonferroni threshold across programs
nom  <- 0.05            # nominal threshold
tb   <- -log10(bonf)    # ~3.08  (signed-score cutoff for Bonferroni)
tn   <- -log10(nom)     # ~1.30  (signed-score cutoff for nominal p<0.05)

sig_bonf <- df$MEANgamma_top100_shet_adjusted_P < bonf | df$P_withShet < bonf
n_bonf   <- sum(sig_bonf)
df$bonf_sig <- sig_bonf
if (n_bonf > 0) {
  ann_lab <- sprintf("Bonferroni-sig: P%s", paste(sort(df$Program[sig_bonf]), collapse = ", P"))
  ann_col <- "#1a7a34"
} else {
  ann_lab <- "No program reached Bonferroni significance"
  ann_col <- "#c0392b"
}

## ---- label programs significant at nominal p<0.05 in EITHER test (Bonferroni fallback) ----
df$LABEL <- as.character(df$Program)
df$LABEL[df$MEANgamma_top100_shet_adjusted_P > nom & df$P_withShet > nom] <- ""

## ---- colour by nominal-p<0.05 category (paper's four-way scheme) ----
df$COLOR <- ifelse(df$MEANgamma_top100_shet_adjusted_P < nom & df$P_withShet <  nom, "both_enriched",
             ifelse(df$MEANgamma_top100_shet_adjusted_P < nom & df$P_withShet >= nom, "program_enriched",
             ifelse(df$MEANgamma_top100_shet_adjusted_P >= nom & df$P_withShet <  nom, "regulator_enriched", "other")))
df$COLOR <- factor(df$COLOR, levels = c("other", "program_enriched", "regulator_enriched", "both_enriched"))

library(ggplot2)
library(ggrepel)

lim <- max(tb, max(abs(df$Program_score)), max(abs(df$regulator_score))) * 1.1

g <- ggplot(df, aes(x = Program_score, y = regulator_score, label = LABEL, color = COLOR)) +
  theme_classic(base_size = 20) +
  # origin
  geom_vline(xintercept = 0, color = "grey85") +
  geom_hline(yintercept = 0, color = "grey85") +
  # nominal p<0.05 cutoffs (grey dashed)
  geom_vline(xintercept = c(-tn, tn), linetype = "dashed", color = "grey55") +
  geom_hline(yintercept = c(-tn, tn), linetype = "dashed", color = "grey55") +
  # Bonferroni cutoffs (red dotted) -- no point crosses these
  geom_vline(xintercept = c(-tb, tb), linetype = "dotted", color = "#c0392b", linewidth = 0.8) +
  geom_hline(yintercept = c(-tb, tb), linetype = "dotted", color = "#c0392b", linewidth = 0.8) +
  geom_point(size = 4.5) +
  # black ring around any Bonferroni-significant program
  geom_point(data = df[df$bonf_sig, , drop = FALSE], shape = 21, size = 8, stroke = 1.4, fill = NA, color = "black") +
  geom_text_repel(size = 6, max.overlaps = Inf, box.padding = 0.5, min.segment.length = 0, segment.color = "grey50") +
  scale_color_manual(
    values = c("other" = "grey70", "program_enriched" = "#FEA601",
               "regulator_enriched" = "#4783B5", "both_enriched" = "#34a853"),
    labels = c("other" = "n.s. (p >= 0.05)", "program_enriched" = "program p<0.05",
               "regulator_enriched" = "regulator p<0.05", "both_enriched" = "both p<0.05"),
    drop = FALSE) +
  annotate("text", x = -lim * 0.98, y = lim * 0.98, hjust = 0, vjust = 1, size = 4.6, fontface = "bold",
           color = ann_col, label = ann_lab) +
  annotate("text", x = tb, y = -lim, angle = 90, hjust = 0, vjust = 1.2, size = 4, color = "#c0392b",
           label = paste0("Bonferroni  0.05/", N)) +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  guides(color = guide_legend(override.aes = aes(label = ""))) +
  theme(legend.title = element_blank(), legend.position = "right",
        plot.title = element_text(size = 19), plot.subtitle = element_text(size = 12, color = "grey30")) +
  labs(
    x = "Program burden effect, signed -log10(P)",
    y = "Regulator-burden correlation, signed -log10(P)",
    title = sprintf("Fig 4C  |  %s  |  %s", LABEL, trait),
    subtitle = sprintf(
      "%s - %d cNMF programs vs %s LoF gamma\nBonferroni 0.05/%d = %.1e (red dotted): %d significant. Nominal p<0.05 (grey dashed): coloured & labelled.",
      LABEL, N, trait, N, bonf, n_bonf),
    caption = "signed -log10(P); black ring = Bonferroni-significant (p < 0.05/60)")

ggsave(plot = g, file.path(OUT, sprintf("Fig4C_%s_%s.pdf", label_slug, trait)), width = 12.5, height = 9)
ggsave(plot = g, file.path(OUT, sprintf("Fig4C_%s_%s.png", label_slug, trait)), width = 12.5, height = 9, dpi = 150)

cat(sprintf("N=%d  Bonferroni=%.3e (-log10=%.2f)  nominal -log10=%.2f\n", N, bonf, tb, tn))
cat(sprintf("Bonferroni-significant programs: %d\n", n_bonf))
cat(sprintf("nominal p<0.05 labelled programs: %s\n",
            paste(sort(df$Program[df$LABEL != ""]), collapse = ", ")))
cat("saved: ", file.path(OUT, sprintf("Fig4C_%s_%s.{pdf,png}", label_slug, trait)), "\n")
