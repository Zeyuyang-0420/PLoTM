#!/usr/bin/env Rscript
# fig5a_regulator_model.R <CT> <TRAIT_FILE> <K> [NTOP=10] [FIGDIR] [SELECTED_txt]
#
# Paper-style Fig 5a for ONE (condition, trait): the regulator -> program -> trait causal chain for
# the RN best-subset regulator programs of the multiple regression
#     gamma_x = a + sum_P w_P beta_{x->P} + lambda*shet .
# Left column = each regulator program's top-NTOP trans-regulators (perturbed genes, BH<0.05, by
# |beta_x(P)|); middle = the regulator programs labelled with w_P; right = the trait.
#   gene->program : sign(beta) up/down-regulates ; width ~ |beta|
#   program->trait: sign(w_P) increases/decreases ; width ~ |w_P|
#   gene node fill = its own LoF gamma (>+thr red / <-thr blue / else grey);  label colour = sign(gamma)
#
# Reads the regulator programs from the SELECTED_txt written by fig5_step1_permutation.R
# (SelectedPrograms_<TAG>.txt) so it depicts exactly the model the permutation test evaluated.

suppressMessages(library(ggplot2))
a <- commandArgs(trailingOnly = TRUE)
CT <- a[1]; TRAIT_FILE <- a[2]; K <- as.numeric(a[3])
NTOP  <- if (length(a) >= 4 && nzchar(a[4])) as.numeric(a[4]) else 10
FIGDIR<- if (length(a) >= 5 && nzchar(a[5])) a[5] else "figures"
SELF  <- if (length(a) >= 6 && nzchar(a[6])) a[6] else NA
LOF_THR <- as.numeric(Sys.getenv("FIG5A_GAMMA_THR", "0.03"))
TRAIT <- sub("\\.per_gene_estimates\\.tsv$", "", TRAIT_FILE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

RESULTS <- Sys.getenv("RESULTS_ROOT", "results")
DATA    <- Sys.getenv("DATA_ROOT", "data")
regdir  <- file.path(RESULTS, "cNMF_regulation", CT)
if (!dir.exists(regdir)) regdir <- file.path(DATA, "Perturbseq", "cNMF_regulation", CT)
cnmfdir <- file.path(RESULTS, "cNMF", CT, "test1")
if (!file.exists(file.path(cnmfdir, sprintf("test1.gene_spectra_score.k_%d.dt_0_4.txt", K))))
  cnmfdir <- file.path(DATA, "Perturbseq", "cNMF", CT, "test1")
ANN <- Sys.getenv("PROGRAM_ANNOTATION", "")     # optional program_annotation_named.tsv (Program,name)

corr <- read.table(file.path(DATA, "gencode_v41_gname_gid_ALL_sorted_onlyID")); corr <- corr[!duplicated(corr[,1]),]; rownames(corr) <- corr[,1]
shet <- read.table(file.path(DATA, "shet_10bins.txt"), header = TRUE)
am <- setNames(paste0("P", 1:K), paste0("P", 1:K))
if (nzchar(ANN) && file.exists(ANN)) { an <- read.table(ANN, header = TRUE, sep = "\t", quote = ""); am <- setNames(an$name, paste0("P", an$Program)) }

## regulator programs = RN best-subset from step1
sel_f <- if (!is.na(SELF)) SELF else file.path(RESULTS, "fig5", CT, "P200", paste0("SelectedPrograms_", TRAIT, "_program5_regulator3_LOF0.03.txt"))
sel <- read.table(sel_f, header = TRUE, sep = "\t")
REG <- strsplit(sel$Regulator_selected[1], ",")[[1]]
cat(sprintf("[fig5a-reg] %s | %s : regulators %s\n", CT, TRAIT, paste(REG, collapse = ",")))

Bb <- NULL; Pp <- NULL
for (p in REG) { i <- as.numeric(sub("^P", "", p))
  x <- read.table(file.path(regdir, paste0("K", K, "_program", i, "_perturb_effects.txt")), header = TRUE)
  b <- x[, c(1,2)]; q <- x[, c(1,3)]; colnames(b) <- c("GENE", p); colnames(q) <- c("GENE", p)
  Bb <- if (is.null(Bb)) b else merge(Bb, b, by = "GENE"); Pp <- if (is.null(Pp)) q else merge(Pp, q, by = "GENE") }
rownames(Bb) <- Bb$GENE; rownames(Pp) <- Pp$GENE

L <- read.table(file.path(DATA, "LoF", "GeneBayes_posterior", TRAIT_FILE), sep = "\t", quote = "", header = TRUE)
L <- data.frame(GENE = corr[L$ensg, 2], L); L <- L[!is.na(L$GENE), ]
gamma <- setNames(L$post_mean, L$GENE)

## program -> trait weights (the RN-program scaled regression, exactly score_model's fit1)
d <- merge(L, Bb, by = "GENE"); d <- merge(d, shet, by = "ensg"); dd <- d[, c("post_mean", REG, "shet")]
for (j in 1:ncol(dd)) { v <- dd[,j]; v[is.infinite(v)] <- max(v[!is.infinite(v)]); dd[,j] <- scale(v) }
co <- summary(lm(post_mean ~ ., data = dd))$coefficients
wP <- co[REG, "Estimate"]; wP_p <- co[REG, "Pr(>|t|)"]

edges <- do.call(rbind, lapply(REG, function(p) {
  padj <- p.adjust(Pp[, p], method = "BH"); sig <- rownames(Pp)[padj < 0.05]
  bs <- Bb[sig, p]; g <- head(sig[order(-abs(bs))], NTOP)
  data.frame(GENE = g, PROGRAM = p, beta = Bb[g, p], stringsAsFactors = FALSE) }))
genes <- sort(unique(edges$GENE))

gy <- setNames(seq_along(genes), genes)
py <- setNames(seq(quantile(gy, .12), quantile(gy, .88), length.out = length(REG)), REG)
ty <- mean(range(gy)); relw <- function(v) pmax(abs(v)/max(abs(v)), 0.12)
ed <- data.frame(x=0, y=gy[edges$GENE], xend=1, yend=py[edges$PROGRAM],
                 dir=ifelse(edges$beta>0,"up-regulates","down-regulates"), w=0.3+1.9*relw(edges$beta))
pt <- data.frame(x=1, y=py[REG], xend=2, yend=ty,
                 dir=ifelse(wP>0,"increases trait gamma","decreases trait gamma"), w=0.6+2.4*relw(wP))
gv <- gamma[genes]; gv[is.na(gv)] <- 0
gcls <- ifelse(gv>LOF_THR,"gamma > +thr", ifelse(gv< -LOF_THR,"gamma < -thr","|gamma| small"))
lcol <- ifelse(gv>0,"#c0392b", ifelse(gv<0,"#2c6fbb","grey30"))
gn <- data.frame(x=0, y=gy[genes], lab=genes, cls=gcls, lcol=lcol, stringsAsFactors=FALSE)
wrap <- function(s,w=22){ss<-strwrap(s,w); paste(ss[1:min(2,length(ss))],collapse="\n")}
pn <- data.frame(x=1, y=py[REG], lab=vapply(REG, function(p) paste0(p,"\n",wrap(am[[p]])), character(1)))
pw <- data.frame(x=1.5, y=(py[REG]+ty)/2, lab=sprintf("w=%+.3f\nP=%.2g", wP, wP_p))

g <- ggplot() + theme_void(base_size = 14) +
  geom_segment(data=ed, aes(x,y,xend=xend,yend=yend,color=dir,linewidth=w), alpha=0.6) +
  geom_segment(data=pt, aes(x,y,xend=xend,yend=yend,color=dir,linewidth=w), alpha=0.9,
               arrow=arrow(length=unit(0.18,"cm"), type="closed")) +
  geom_point(data=gn, aes(x,y,fill=cls), shape=21, size=3.4, color="grey25") +
  geom_text(data=gn, aes(x-0.03,y,label=lab), color=gn$lcol, hjust=1, size=2.9) +
  geom_label(data=pn, aes(x,y,label=lab), size=3.4, fill="grey92", lineheight=0.9, label.padding=unit(0.4,"lines")) +
  geom_text(data=pw, aes(x,y,label=lab), size=3.0, color="grey25", lineheight=0.9) +
  geom_label(data=data.frame(x=2,y=ty), aes(x,y), label=TRAIT, size=4.4, fill="#ffe3e3", fontface="bold") +
  scale_linewidth_identity() +
  scale_color_manual(name=NULL, values=c("up-regulates"="#c0392b","down-regulates"="#2c6fbb",
    "increases trait gamma"="#c0392b","decreases trait gamma"="#2c6fbb")) +
  scale_fill_manual(name="regulator LoF", values=c("gamma > +thr"="#e8746a","gamma < -thr"="#6a9fe8","|gamma| small"="grey85")) +
  guides(color=guide_legend(override.aes=list(linewidth=1.4), nrow=2, order=2),
         fill=guide_legend(override.aes=list(size=4), nrow=2, order=1)) +
  coord_cartesian(xlim=c(-0.55,2.25), clip="off") +
  labs(title=sprintf("Fig 5a (regulator model) | %s | %s", CT, TRAIT),
       subtitle=sprintf(paste0("regulators %s -> trait. Left = top-%d trans-regulators (BH<0.05, |beta|). ",
         "gene label colour = sign(gamma).\nEdge width: gene->program |beta| ; program->trait |w_P|. program->trait P: %s"),
         paste(REG, collapse=", "), NTOP, paste(sprintf("%s %.2g", REG, wP_p), collapse=" ; "))) +
  theme(legend.position="bottom", legend.box="horizontal", legend.text=element_text(size=9),
        plot.title=element_text(size=14,face="bold"), plot.subtitle=element_text(size=8.6,color="grey30"),
        plot.background=element_rect(fill="white",colour=NA), plot.margin=margin(10,14,8,26))
h <- max(7, length(genes)*0.34)
out <- file.path(FIGDIR, sprintf("Fig5a_regulatorModel_%s_%s", TRAIT, CT))
ggsave(paste0(out,".png"), g, width=13, height=h, dpi=150, limitsize=FALSE, bg="white")
ggsave(paste0(out,".pdf"), g, width=13, height=h, limitsize=FALSE, bg="white")
write.table(edges, paste0(out,"_edges.tsv"), row.names=FALSE, sep="\t", quote=FALSE)
cat(sprintf("[fig5a-reg] %d genes -> %s.{png,pdf}\n", length(genes), out))
