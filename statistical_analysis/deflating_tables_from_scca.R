#### read in data ####
snp <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twins240_noissHweLD_SNPresiduals.csv", header=T, sep=','), row.names=1)
rownames(snp) <- snp[,1]
snp[,1] <- NULL
gene <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twins240_noissHweLD_GENEresiduals.csv", header=T, sep=','), row.names=1)
rownames(gene) <- gene[,1]
gene[,1] <- NULL

alpha <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_alpha_results4.5.csv")
beta <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_beta_results4.5.csv")

#### ensure proper permutation of people ####
sum(colnames(snp)!=alpha[,1])
sum(colnames(gene)!=beta[,1])

#### deflate ####
snp.mat <- as.matrix(snp)
u <- snp.mat %*% alpha[,2]
snp.deflated <- snp.mat - u %*% (t(u) %*% snp.mat)/(sum(u*u))
write.csv(snp.deflated, "/workdir/users/fnn3/twins_uk/scca/twins240_deflatedSNPs_scca4.5.csv",  quote=F)

gene.mat <- as.matrix(gene)
v <- gene.mat %*% beta[,2]
gene.deflated <- gene.mat - v %*% (t(v) %*% gene.mat)/(sum(v*v))
write.csv(gene.deflated, "/workdir/users/fnn3/twins_uk/scca/twins240_deflatedGenes_scca4.5.csv", quote=F) 