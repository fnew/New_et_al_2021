#### read in data ####
library(data.table)
snp <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/snp_residuals_mod.csv", header=T, sep=','), row.names=1)
rownames(snp) <- snp[,1]
snp[,1] <- NULL
gene <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/abundance_residuals_mod.csv", header=T, sep='\t'), row.names=1)
rownames(gene) <- gene[,1]
gene[,1] <- NULL

alpha <- read.csv("/workdir/users/fnn3/twins_uk/scca/gene_results/twinsUK_alpha_results.csv")
beta <- read.table("/workdir/users/fnn3/twins_uk/scca/gene_results/twinsUK_beta_results.csv", sep=",", header=T)

#### ensure proper permutation of people ####
sum(colnames(snp)!=alpha[,1])
sum(colnames(gene)!=beta[,1])

#### deflate ####
snp.mat <- as.matrix(snp)
u <- snp.mat %*% alpha[,2]
snp.deflated <- snp.mat - u %*% (t(u) %*% snp.mat)/(sum(u*u))
snp.deflated <- as.data.frame(snp.deflated)
write.csv(snp.deflated, "/workdir/users/fnn3/twins_uk/scca/gene_results/twins240_deflatedSNPs_scca.csv",  quote=F, row.names=F)

gene.mat <- as.matrix(gene)
v <- gene.mat %*% beta[,2]
gene.deflated <- gene.mat - v %*% (t(v) %*% gene.mat)/(sum(v*v))
gene.deflated <- as.data.frame(gene.deflated)
write.table(gene.deflated, "/workdir/users/fnn3/twins_uk/scca/gene_results/twins240_deflatedGenes_scca.csv", quote=F, sep=",", row.names=F) 


##########Species
snp <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/snp_residuals_mod.csv", header=T, sep=','), row.names=1)
rownames(snp) <- snp[,1]
snp[,1] <- NULL
spec <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/species/species_abundance_residuals_v2.csv", header=T, sep='\t'), row.names=1)
rownames(spec) <- spec[,1]
spec[,1] <- NULL

alphaS <- read.csv("/workdir/users/fnn3/twins_uk/scca/species/twinsUK_species_alpha_results_v2.csv")
betaS <- read.csv("/workdir/users/fnn3/twins_uk/scca/species/twinsUK_species_beta_results_v2.csv", header=T)

#### ensure proper permutation of people ####
sum(colnames(snp)!=alphaS[,1])
sum(colnames(spec)!=betaS[,1])

#### deflate ####
snp.mat <- as.matrix(snp)
u <- snp.mat %*% alphaS[,2]
snp.deflatedS <- snp.mat - u %*% (t(u) %*% snp.mat)/(sum(u*u))
snp.deflatedS <- as.data.frame(snp.deflatedS)
write.csv(snp.deflatedS, "/workdir/users/fnn3/twins_uk/scca/species/twins240_deflatedSNPs_species_v2.csv",  quote=F, row.names=F) #Make sure it's the species one

spec.mat <- as.matrix(spec)
v <- spec.mat %*% betaS[,2]
spec.deflated <- spec.mat - v %*% (t(v) %*% spec.mat)/(sum(v*v))
spec.deflated <- as.data.frame(spec.deflated)
write.table(spec.deflated, "/workdir/users/fnn3/twins_uk/scca/species/twins240_deflatedSpecies_v2.csv", quote=F, sep="\t", row.names=F) 

