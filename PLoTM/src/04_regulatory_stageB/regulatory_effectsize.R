args   <- commandArgs(trailingOnly = T)

options("scipen"=10) ##disable 100000 to 1e6

CELL<-as.character( args[1] ) ##
i <-as.numeric( args[2] ) ##1-60
## BACKWARD-COMPATIBLE ADDITIONS (default to the pure-Perturb-seq behaviour):
##   args[3]=K            cNMF program count (default 60). The disease-single-cell basis uses K=30.
##   args[4]=min_profiles per-gene profile gate (default 4). Unchanged for pure Perturb-seq.
K <- if (length(args) >= 3 && nzchar(args[3])) as.numeric(args[3]) else 60
MINP <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else 4
i=i + 1 ##adjust for the first column

metadata=read.csv(paste0("data/Perturbseq/metadata/", CELL, "_metadata.csv"), row.names = 1)

library(data.table)
count<-fread(paste0("data/Perturbseq/cNMF/", CELL, "/test1/test1.usages.k_", K, ".dt_0_4.consensus.txt"), header=T, data.table=F)

## ESTIMABILITY GUARD (no-op for pure Perturb-seq; needed for projected disease-basis usages):
## a program whose usage has zero variance among the non-targeting controls admits no perturb-vs-
## control contrast -> lm() dies with "contrasts can be applied only to factors with 2 or more
## levels". Emit an empty beta file (so downstream merges stay aligned) and stop.
ctrl0 <- row.names(metadata)[metadata$gene == "non-targeting"]
if (sd(count[is.element(count$cell_barcode, ctrl0), i], na.rm=TRUE) == 0) {
  dir.create(paste0("data/Perturbseq/cNMF_regulation/", CELL), showWarnings=FALSE, recursive=TRUE)
  write.table(data.frame(GENE=character(0), lm_es=numeric(0), lm_p=numeric(0)),
              paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K", K, "_program", colnames(count)[i], "_perturb_effects.txt"),
              row.names=FALSE, sep="\t", quote=FALSE)
  cat("SKIP program", colnames(count)[i], "- zero usage variance among NTC controls (not estimable)\n")
  quit(save="no", status=0)
}

GENES<-unique(sort(metadata$gene))
GENES<-GENES[!is.element(GENES, "non-targeting")]

library(dplyr)

control<- row.names( metadata %>% filter(gene %in% "non-targeting") )

N=ncol(count)

summary<-data.frame()

B<- count[is.element(count$cell_barcode, control),c(1,i)]
dfB<-cbind(B, metadata[B$cell_barcode, c("gem_group", "n_genes", "mitopercent")])
dfB<-data.frame(dfB, GROUP="control")

A<- count[,c(1,i)]
tmp<-intersect(unique(A$cell_barcode), row.names(metadata))
row.names(A)<-A$cell_barcode
A<-A[tmp,]
A<-cbind(A, metadata[ tmp, c("gem_group", "n_genes", "mitopercent")])



for(GENE in GENES) {
target<- row.names( metadata %>% filter(gene %in% GENE ) )
## PSEUDOBULK MODALITY ADAPTATION (documented, non-scientific): the paper's single-cell
## gate is ">10 cells" per perturbation. Here each observation is a PSEUDOBULK profile
## (one per donor x guide within Stim48hr, which exists only in run R2), so a gene has a
## median of ~7 profiles and >10 leaves only 1 gene testable. Loosen the gate to >4 profiles
## so the same lm (vs 3,441 non-targeting controls) can be fit. The model, covariates and
## effect-size definition are UNCHANGED; only the applicability threshold is adapted to the
## pseudobulk replicate structure.
if(length(target)>MINP){

dfA<-A[is.element(A$cell_barcode, target),]
dfA<-data.frame(dfA, GROUP="perturb")

if(length(which(is.element(dfB$gem_group, unique(dfA$gem_group))))>0){

df<-rbind(dfA,dfB)
df$gem_group <-factor(as.character(df$gem_group))

colnames(df)[2]<-"exp"
## safety net (never fires in pure Perturb-seq; the NTC guard above already ensures variance):
## skip a gene whose subset usage is constant or leaves a single gem_group level.
if (sd(df$exp, na.rm=TRUE) == 0 || length(unique(df$gem_group)) < 2) next
df$exp<-scale(df$exp)

m.lm <- lm(exp ~ GROUP  + gem_group + n_genes + mitopercent, data=df)
lm_es<-summary(m.lm)$coefficients[2,1]
lm_p<-summary(m.lm)$coefficients[2,4]
hoge<-data.frame(GENE=GENE,  lm_es=lm_es, lm_p=lm_p)
summary<-rbind(summary, hoge)

}}}

dir.create(paste0("data/Perturbseq/cNMF_regulation/", CELL), showWarnings=F)
write.table(summary, paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K", K, "_program", colnames(count)[i],  "_perturb_effects.txt"), row.names=F, sep="\t", quote=F, append=F)
