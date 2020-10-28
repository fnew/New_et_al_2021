setwd("/workdir/users/fnn3/scripts/snp_metagenomes/statistical_analysis")

#Load libraries
library("parallel")
library("doParallel")
library("optparse")
library("data.table")
library("lme4")


#Load functions
#sourcce("./get_scca_functions.R")
source("./refactored_cca_functions.r")

#Set command line options
option_list=list(
	make_option(c("-s", "--snp"), action="store", type="character", default=NULL, help="Input the SNP residuals", metavar="character"),
	make_option(c("-g", "--gene"), action="store", type="character", default=NULL, help="Input the gene abundance residuals", metavar="character"),
	make_option(c("-a", "--alpha"), action="store", type="character", default=NULL, help="Output ALPHA results from the model", metavar="character"),
	make_option(c("-b", "--beta"), action="store", type="character", default=NULL, help="Output BETA results from the model", metavar="character"),
	make_option(c("-c", "--cv"), action="store", type="character", default=NULL, help="Output CV results from the model for component 1", metavar="character"),
	make_option(c("-m", "--cv2"), action="store", type="character", default=NULL, help="Output CV results from the model for component 2", metavar="character"),
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
p = ncol(snp.residuals)/2
group <- rep(1:p, each = 2)


obj <- scca_front(snp.residuals, ab.residuals, 
            penalization_x = "glasso", penalization_y = "enet",
            num_components=2,
            group_x = group, 
            cross_validate=TRUE, parallel_CV = FALSE,
            nlambda=12)


#Output: we want ALPHA, BETA, CV_RESULTS
fwrite(data.frame(obj$alpha), opt$alpha, row.names=T,quote=F)
fwrite(data.frame(obj$beta), opt$beta, row.names=T,quote=F)
fwrite(data.frame(obj$all_other[[1]][[1]]$abs_cors), opt$cv, quote=F)
fwrite(data.frame(obj$all_other[[2]][[1]]$abs_cors), opt$cv2, quote=F)
