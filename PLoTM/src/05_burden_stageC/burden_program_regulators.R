args   <- commandArgs(trailingOnly = T)
options("scipen"=10) ##disable 100000 to 1e6
FILE<-as.character( args[1] ) ##e.g., Backman_2021_86.per_gene_estimates.tsv
K<-as.numeric( args[2] ) ##60
## REPRODUCTION FIX 1 (documented, non-scientific): the original script references
## `CELL` (for the essential-cNMF file paths) but never reads it from args, so it
## errors with "object 'CELL' not found". Read it as the 3rd argument, matching the
## CT name used in Stage A/B (e.g. K562essential_subset). No change to any statistic.
CELL<-as.character( args[3] ) ##e.g., K562essential_subset
## HARNESS PARALLELIZATION (non-scientific): process a SINGLE program per invocation so the
## 100k-permutation program-burden test can run across programs on many cores. args[4] is the
## 1-based program index. All statistics (cor.test, lm, the set.seed()-based permutation) are
## byte-identical to the original; only the loop bound and output paths change. Per-program
## outputs are concatenated back into the standard regulators_/programs_enrichment files after.
PROG<-as.numeric( args[4] )

##data loading
## REPRODUCTION FIX 2 (documented, non-scientific): the essential-cNMF consensus in
## Stage A runs at --local-density-threshold 0.4 -> files are written as dt_0_4, but
## the original path hardcodes dt_0_5 (the genome-wide setting). Point it at the
## dt_0_4 file that Stage A actually produces. Same file, same content otherwise.
GEP<-read.table(paste0("data/Perturbseq/cNMF/", CELL, "/test1/test1.gene_spectra_score.k_", K, ".dt_0_4.txt"), header=T, stringsAsFactor=F)
GEP<-t(GEP)
colnames(GEP)<-paste0("P", c(1:K))
corresp<-read.table("data/gencode_v41_gname_gid_ALL_sorted_onlyID", header=F, stringsAsFactor=F)
corresp<-corresp[!duplicated(corresp[,1]),]
row.names(corresp)<-corresp[,1]
GEP<-data.frame(GENE=corresp[row.names(GEP), 2], GEP)
## REPRODUCTION FIX 3 (documented, non-scientific; same as Figure4/1): the regulator loop
## maps each perturbed gene SYMBOL to its ENSG to merge with the ENSG-keyed LoF table. The
## shipped `corresp[as.character(tmp$GENE),1]` indexes by symbol while row.names are ENSG ->
## all NA -> empty merge -> cor.test("not enough finite observations"). Add a symbol-keyed
## lookup so the SAME symbol->ENSG mapping resolves. Verified to reproduce the genome-wide
## regulators_enrichment numbers exactly. No statistic changed.
corresp_bySym<-corresp[!duplicated(corresp[,2]),]
row.names(corresp_bySym)<-corresp_bySym[,2]

## ESTIMABILITY (no-op for pure Perturb-seq): with the disease-single-cell basis some programs are
## non-estimable (empty beta file). Restrict the beta merge to the ALIVE programs listed in the file
## named by env ALIVE_PROGRAMS; absent -> all 1:K, so pure-Perturb-seq behaviour is identical.
.alivef <- Sys.getenv("ALIVE_PROGRAMS")
if (!nzchar(.alivef)) for (cand in c(file.path(Sys.getenv("RESULTS_ROOT",""), "cNMF", CELL, "alive_programs.txt"),
                                     file.path("data/Perturbseq/cNMF", CELL, "alive_programs.txt")))
  if (file.exists(cand)) { .alivef <- cand; break }
PROGS <- if (nzchar(.alivef) && file.exists(.alivef)) as.numeric(readLines(.alivef)) else 1:K
PROG_ALIVE <- PROG %in% PROGS

GEP_reg_beta<-data.frame()
.first<-TRUE
for(i in PROGS){
tmp<-read.table(paste0("data/Perturbseq/cNMF_regulation/", CELL, "/K", K, "_program", i, "_perturb_effects.txt"), header=T)
tmp1<-tmp[,c(1,2)]

colnames(tmp1)<-c("GENE", paste0("P", i))

if(.first){
GEP_reg_beta<-tmp1; .first<-FALSE
} else {
GEP_reg_beta<-merge(GEP_reg_beta, tmp1, by="GENE")
}
}

LOF<-read.table(paste0("data/LoF/GeneBayes_posterior/", FILE), sep="\t", quote="", header=T, stringsAsFactor=F)

LOF$post_mean[is.element(LOF$post_mean, "Inf")]<-max(LOF$post_mean[!is.infinite(LOF$post_mean)])
LOF$post_mean[is.element(LOF$post_mean, "-Inf")]<-min(LOF$post_mean[!is.infinite(LOF$post_mean)])

LOF<-data.frame(LOF, gene=corresp[as.character(LOF$ensg),2])

GEP2<-data.frame(ensg=row.names(GEP), GEP)
df<-merge(GEP2, LOF, by="ensg")


shet<-read.table("data/shet_10bins.txt", header=T, stringsAsFactor=F)
shet<-shet[is.element(shet$ensg, LOF$ensg),]

summary_reg<-data.frame()
summary_pro<-data.frame()

for(Program in c(PROG)){

if (PROG_ALIVE) {
tmp<-GEP_reg_beta[,c("GENE", paste0("P", Program))]
colnames(tmp)<-c("GENE", "perturb_beta")
tmp<-data.frame(tmp, ensg=corresp_bySym[ as.character(tmp$GENE), 1])  ##FIX 3: symbol->ENSG
df<-merge(tmp, LOF, by="ensg")

P_pearson_all<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$p.value
R_pearson_all<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$estimate
R_pearson_all_CIlower<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[1]
R_pearson_all_CIupper<-cor.test(df$perturb_beta, df$post_mean, method="pearson")$conf.int[2]

##with shet regression
df<-merge(df, shet, by="ensg")
df$post_mean<-scale(df$post_mean)
df$perturb_beta<-scale(df$perturb_beta)

fit<- lm(post_mean~perturb_beta + shet, data=df)
P_withShet<-summary(fit)$coefficients[2,4]
beta_withShet<-summary(fit)$coefficients[2,1]
betaSE_withShet<-summary(fit)$coefficients[2,2]
} else {
## non-estimable program (no beta_x(P)): emit NA regulator row; program-burden below still runs.
P_pearson_all<-NA; R_pearson_all<-NA; R_pearson_all_CIlower<-NA; R_pearson_all_CIupper<-NA
P_withShet<-NA; beta_withShet<-NA; betaSE_withShet<-NA
}

hoge<-data.frame(FILE=FILE, Program=Program,
    P_pearson_all=P_pearson_all, R_pearson_all=R_pearson_all,
    R_pearson_all_CIlower=R_pearson_all_CIlower, R_pearson_all_CIupper=R_pearson_all_CIupper,
    P_withShet=P_withShet, beta_withShet=beta_withShet, betaSE_withShet=betaSE_withShet)
summary_reg<-rbind(summary_reg, hoge)

###
##from here, program enrichment

##1. mean(gamma) for top 100 program genes

tmp<-GEP[,c("GENE",paste0("P", Program))]
tmp<-tmp[order(tmp[,2], decreasing=T),]
colnames(tmp)[2]<-"GEPscore"
tmp<-data.frame(tmp, ensg=row.names(tmp))
df<-merge(tmp, LOF, by="ensg")
df2<-df[order(df$GEPscore, decreasing=T),][1:100,]

MEANgamma_top100=mean(df2[,"post_mean"])

 
shet_tmp<-shet[is.element(shet$ensg, df2$ensg),]
shet_tmp_all<-shet[is.element(shet$ensg, df$ensg),]

A<-table(shet_tmp$shet_BIN)

random<-c()
for(I in 1:100000){
set.seed(I)
genes<-c()
for(p in 1:length(A)){
genes<-c(genes, sample(shet_tmp_all$ensg[is.element(shet_tmp_all$shet_BIN, names(A)[p])], A[p]))
}
random<-c(random, mean(df[is.element(df$ensg, genes),"post_mean"], na.omit=T))
}

random_mean=mean(random)


P1=(rank(c(MEANgamma_top100, random))[1]/length(random)) * 2
P2=(( length(random) + 2 - (rank( c(MEANgamma_top100, random))[1]) ) /length(random)) * 2
program_P=min(P1,P2)


hoge<-data.frame(FILE=FILE, Program=Program, 
MEANgamma_top100=MEANgamma_top100,
shet_adjusted_random_mean=random_mean,
MEANgamma_top100_shet_adjusted_P=program_P
)

summary_pro<-rbind(summary_pro, hoge)

}

dir.create(paste0("data/Perturbseq/trait_association/", CELL, "/ProgramLevel/tmp_parallel"), showWarnings=F, recursive=T)
write.table(summary_reg, paste0("data/Perturbseq/trait_association/", CELL, "/ProgramLevel/tmp_parallel/regulators_enrichment_K", K, "_", FILE, ".p", PROG), row.names=F, sep="\t", quote=F)
write.table(summary_pro, paste0("data/Perturbseq/trait_association/", CELL, "/ProgramLevel/tmp_parallel/programs_enrichment_K", K, "_", FILE, ".p", PROG), row.names=F, sep="\t", quote=F)
