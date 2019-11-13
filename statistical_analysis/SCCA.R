setwd("/workdir/users/fnn3/scripts/snp_metagenomes/statistical_analysis")

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
	make_option(c("-m", "--mac"), action="store", type="character", default=NULL, help="Output Mean Abs Corrs results from the model", metavar="character"),
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


#obj <- runCCA(snp.residuals, ab.residuals, 
#              penalization_x = "glasso", penalization_y = "enet", 
#              group_x = group,
#              parallel_CV = F,
#	      nlambda=10)  
              #nr_cores = 16) #Runs fine on one core, not going to use this 

#runCCA is to get the tuning parameters, sCCA takes the tuning parameters that we pick to run

# these are the tuning parameter pairs from the second run
#nnx <- c(0.00963166838580798, 0.007702851, 0.00577403350158676, 0.004617738, 0.0034614421450152, 0.002768261, 0.00207508004932683, 0.00124397780772249, 0.000745745103475914)
#nny <- c(0.0112327319448012, 0.008983289, 0.00673384588896276, 0.00538534, 0.00403683455450812, 0.003228427, 0.00242001873657089, 0.00145076311805099, 0.00086970964021517, 0.000521377231659547)

### this is the tuning parameter pair from the third run (small, minor run)
#nnx <- 0.007702851
#nny <- 0.01123273

# this is the tuning parameter pairs from the FOURTH run (nny is new, nnx is as from 2nd run)
nny <- c(0.0100452415261994, 0.008983289, 0.00803360288023392, 0.0112327319448012, 
  0.008983289, 0.00673384588896276, 0.00538534, 0.00403683455450812, 
  0.003228427, 0.00242001873657089)
nnx <- c(0.00963166838580798, 0.007702851, 0.00577403350158676, 0.004617738, 0.0034614421450152, 0.002768261, 0.00207508004932683, 0.00124397780772249, 0.000745745103475914)


obj <- sCCA(snp.residuals, ab.residuals, 
            penalization_x = "glasso", penalization_y = "enet",
            group_x = group,
            grp_penalty_x = nnx, lasso_penalty_y = nny, 
            cross_validate=TRUE)

#Output: we want ALPHA, BETA, CV_RESULTS
fwrite(data.frame(obj$ALPHA), opt$alpha, row.names=T,quote=F)
fwrite(data.frame(obj$BETA), opt$beta, row.names=T,quote=F)
fwrite(data.frame(obj$CV_results$mean_abs_cors), opt$mac, quote=F)
fwrite(data.frame(obj$CV_results), opt$cv, quote=F)
