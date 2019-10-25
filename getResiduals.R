##########################################################################################################
### This script will extract the residuals from the SNP and abundance tables for use in the SCCA #########
##########################################################################################################

setwd("/workdir/users/fnn3/scripts/twins_uk/scca/")

#Load libraries
library("parallel")
library("doParallel")
library("optparse")
library("data.table")
library("lme4")

#Set command line options######################################################################################
option_list=list(
	make_option(c("-s", "--snp"), action="store", type="character", default=NULL, help="Input the SNP file", metavar="character"),
	make_option(c("-g", "--gene"), action="store", type="character", default=NULL, help="Input the gene abundance table", metavar="character"),
	make_option(c("-m", "--meta"), action="store", type="character", default=NULL, help="Input the metadata table", metavar="character"),
	make_option(c("-r", "--residuals1"), action="store", type="character", default=NULL, help="Output the residuals from the SNP table", metavar="character"),
	make_option(c("-e", "--residuals2"), action="store", type="character", default=NULL, help="Output the residuals from the abundance table", metavar="character"));
opt_parser=OptionParser(option_list=option_list);
opt=parse_args(opt_parser);


#loading data###############################################################################################
X1.a <- as.data.frame(fread(opt$snp, header=T, sep=','), row.names=1)  
rownames(X1.a) <- X1.a[,1]
X1.a[,1] <- NULL

X2.a <- as.data.frame(fread(opt$gene, header=T, sep=','), row.names=1)
rownames(X2.a) <- X2.a[,1]
X2.a[,1] <- NULL

#Import metadata
meta <- read.csv(opt$meta, header=T, sep=',', row.names=1)

#Sort all datasets by rowname
X1.as <- X1.a[ order(row.names(X1.a)), ]
X2.as <- X2.a[ order(row.names(X2.a)), ]
meta.s <- meta[ order(row.names(meta)), ]

#Projecting out the covariates ########################
#Do PCA using only one person per family and the lonely twins
names <- c("NA5172", "NA5369", "NA5937", "NA5530", "NA9465", "NA9457", "NA9419", "NA9229") 
lonely <- (meta.s)[((meta.s$IndividualFamilyID) %in% names),]
pair <- meta.s[!((meta.s$IndividualFamilyID) %in% names),]

toDelete <- seq(1, nrow(pair), 2)
singles <- pair[ toDelete ,]

all.single <- rbind(singles, lonely) #n=124
tokeep <- rownames(all.single)
snp.sing <- subset(X1.as, rownames(X1.as) %in% tokeep)

snp.pca <- prcomp(snp.sing)

#Use pca.predict() to get the PCs for the full SNP data, n=240
pred <- predict(snp.pca, X1.as)
pcs <- pred$x[,c(1:10)]

#New: mixed model from Ben to get the residuals
#Make twinstatus variable
lonely <- names(which(table(meta.s$IndividualFamilyID)==1))
lonely <- sapply(lonely, FUN = function(k) which(meta.s$IndividualFamilyID == k))

twinstatus <- ifelse(meta.s$MZTwin==1, "MZ", "DZ")
twinstatus[lonely] <- "lonely"

#The function to get residuals
get.residuals <- function(y) {
      mod <- lmer(y ~ 1
                    + Age.at.metagenomics.sample
		    + IndividualBMI
		    + pcs
		    + (1|SampleShipmentNum)
		    + (0+dummy(twinstatus,"MZ")|IndividualFamilyID) 
		    + (0+dummy(twinstatus,"DZ")|IndividualFamilyID)
		    , data = meta.s)  
residuals(mod)
}

get.snpresiduals <- function(y) {
      mod <- lmer(y ~ 1
      		    + Age.at.metagenomics.sample
		    + pcs
		    + IndividualBMI
		    + (0+dummy(twinstatus,"MZ")|IndividualFamilyID) 
		    + (0+dummy(twinstatus,"DZ")|IndividualFamilyID)
		    , data = meta.s)  
residuals(mod)
}

#Using the function
snp.residuals <- apply(X1.as, 2, get.snpresiduals)
ab.residuals <- apply(X2.as, 2, get.residuals)

#Output the residuals for future use
fwrite(snp.residuals, opt$residuals1, quote=F, sep=',', row.names=T, col.names=T)
fwrite(ab.residuals, opt$residuals2, quote=F, sep=',', row.names=T, col.names=T)

