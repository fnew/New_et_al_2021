setwd("/workdir/users/fnn3/scripts/twins_uk/scca/")

#Load libraries
library("parallel")
library("doParallel")
library("optparse")
library("data.table")
library("lme4")

#Load functions
source("./get_scca_functions.R")


#Set command line options
option_list=list(
	make_option(c("-s", "--snp"), action="store", type="character", default=NULL, help="Input the SNP residuals", metavar="character"),
	make_option(c("-g", "--gene"), action="store", type="character", default=NULL, help="Input the gene abundance residuals", metavar="character"),
	make_option(c("-a", "--alpha"), action="store", type="character", default=NULL, help="Output ALPHA results from the model", metavar="character"),
	make_option(c("-b", "--beta"), action="store", type="character", default=NULL, help="Output BETA results from the model", metavar="character"),
	make_option(c("-c", "--cv"), action="store", type="character", default=NULL, help="Output CV results from the model", metavar="character"),
	make_option(c("-t", "--thread"), action="store", type="integer", default=1, help="How many threads to use"));
opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);


#Import SNP and abundance residuals
X1.a <- as.data.frame(fread(opt$snp, header=T, sep=','), row.names=1)  
rownames(X1.a) <- X1.a[,1]
X1.a[,1] <- NULL

X2.a <- as.data.frame(fread(opt$gene, header=T, sep=','), row.names=1)
rownames(X2.a) <- X2.a[,1]
X2.a[,1] <- NULL


#Sort all datasets by rowname
snp.residuals <- X1.a[ order(row.names(X1.a)), ]
ab.residuals  <- X2.a[ order(row.names(X2.a)), ]


#Code to get output
group <- rep(seq(from = 1, to = ncol(snp.residuals)/2), 2)


obj <- runCCA(snp.residuals, ab.residuals, 
              penalization_x = "glasso", penalization_y = "enet", 
              group_x = group,
              parallel_CV = F,
	      nlambda=10)  
              #nr_cores = 16) #Runs fine on one core, not going to use this 

#obj <- sCCA(snp.residuals, ab.residuals, lasso_penalty_x = 1/50, lasso_penalty_y = 1/50)

#Output: we want ALPHA, BETA, CV_RESULTS
fwrite(data.frame(obj$ALPHA), opt$alpha, row.names=T,quote=F)
fwrite(data.frame(obj$BETA), opt$beta, row.names=T,quote=F)
fwrite(data.frame(obj$CV_results), opt$cv,row.names=T, quote=F)
