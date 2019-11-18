ndecs <- 4

## read in data matrices and relevant output from runCCA ########################

snp <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twins240_noissHweLD_SNPresiduals.csv", header=T, sep=','), row.names=1)
rownames(snp) <- snp[,1]
snp[,1] <- NULL
gene <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twins240_noissHweLD_GENEresiduals.csv", header=T, sep=','), row.names=1)
rownames(gene) <- gene[,1]
gene[,1] <- NULL

alpha <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_alpha_results.csv")
beta <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_beta_results.csv")
cv <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twinsUK_cv_results.csv", header=T, sep=','), row.names=1)

## inspect the output ########################################################

mat <- cv[, 102:108]
nlambda <- sqrt(nrow(mat)) #may not always work due the coercion that happening during saving

m <- matrix(mat[,7], nrow=nlambda, byrow=FALSE)
rownames(m) <- round(mat[1:nlambda,2], ndecs)
colnames(m) <- round(mat[1+(0:(nlambda-1))*nlambda, 6], ndecs)

# visualize the values
gplots::heatmap.2(m, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

# calculating the max and min
mat[c( which.min(mat[,7]), which.max(mat[,7]) ), ]

# the reported max and min
cv[1,(ncol(cv)-5):ncol(cv)]

## findings #################################################################
# the high area is at the top left of plot, so we need to try: bigger vals of x penalty and of y penalty

# previous penalty params
pen.y <- sort(unique(mat[,2]), dec=TRUE)
pen.x <- sort(unique(mat[,6]), dec=TRUE)

# new ones
dx <- log(pen.x[1]) - log(pen.x[2])
new.pen.x <- c(exp(log(pen.x)[1]+dx*(3:1)), pen.x[1:3]) #3 news ones, 3 old ones
#c(0.00963166838580798, 0.00577403350158676, 0.0034614421450152, 0.00207508004932683, 0.00124397780772249, 0.000745745103475914)

dy <- log(pen.y[1]) - log(pen.y[2])
new.pen.y <- c(exp(log(pen.y)[1]+dy*(3:1)), pen.y[1:4]) #3 new ones, 4 old ones
#c(0.0112327319448012, 0.00673384588896276, 0.00403683455450812, 0.00242001873657089, 0.00145076311805099, 0.00086970964021517, 0.000521377231659547)


####################################################################################################################
### sCCA 2 results
####################################################################################################################

alpha2 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_alpha_results2.csv")
beta2 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_beta_results2.csv")
mac <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_mean_abs_corrs_results2.csv")
cv <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twinsUK_cv_results2.csv", header=T, sep=','), row.names=1)

# make matrix of correlation values for plotting
ndecs <- 4
nlambda_y <- 10; nlambda_x <- 9
m <- matrix(mac[,8], nrow=nlambda_y, byrow=FALSE) #arbitrary test
rownames(m) <- round(mac[1:nlambda_y,3], ndecs)
colnames(m) <- round(mac[1+(0:(nlambda_x-1))*nlambda_y, 7], ndecs)

# visualize the values
gplots::heatmap.2(m, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

# Min and Max
mac[c( which.min(mac[,8]), which.max(mac[,8]) ), ]

# the reported max and min
cv[2,(ncol(cv)-6):ncol(cv)]

############
#easier visualization for 1 se rule
OneSEOfMax <- mac[which.max(mac[,8]), 8] - mac[which.max(mac[,8]), 9]

m.thresholded <- m
m.thresholded[mac$mean.Cor.over.CVs<=OneSEOfMax] <- NA_real_

gplots::heatmap.2(m.thresholded, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

### looking at this plot, the "top left" point is:
lambda_y_1se <- rev(sort(unique(mac[,3])))[1] #1 chosen from looking at plot
#0.01123273
lambda_x_1se <- rev(sort(unique(mac[,7])))[2] #2 chosen from looking at plot
#0.007702851


#############################################################
## Third run of SCCA
#############################################################

alpha3 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_alpha_results3.csv")
#356
beta3 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_beta_results3.csv")
#171

#################################################################
## Trying a few more lambdas (y_lambdas) for the fourth run

k <- 3 #the number of "new" values
d <- log(nny)[1] - log(nny)[2]
nny2 <- c(-(1:k)*(0.5*d)+log(nny)[1], log(nny)[1:(length(nny)-k)])
nny2 <- exp(nny2)

###################################################################
## Fourth run of SCCA 
#####################################################################

alpha4 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_alpha_results4.csv")
beta4 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_beta_results4.csv")
mac4 <- read.csv("/workdir/users/fnn3/twins_uk/scca/twinsUK_mean_abs_corrs_results4.csv")
#cv4 <- as.data.frame(fread("/workdir/users/fnn3/twins_uk/scca/twinsUK_cv_results4.csv", header=T, sep=','), row.names=1)

# make matrix of correlation values for plotting
ndecs <- 4
nlambda_y <- 10; nlambda_x <- 9
m <- matrix(mac4[,7], nrow=nlambda_y, byrow=FALSE) #arbitrary test
rownames(m) <- round(mac4[1:nlambda_y,2], ndecs)
colnames(m) <- round(mac4[1+(0:(nlambda_x-1))*nlambda_y, 6], ndecs)

# visualize the values
gplots::heatmap.2(m, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
# Min and Max
mac4[c( which.min(mac4[,7]), which.max(mac4[,7]) ), ]


#easier visualization for 1 se rule
OneSEOfMax <- mac4[which.max(mac4[,7]), 7] - mac4[which.max(mac4[,7]), 8]

m.thresholded <- m
m.thresholded[mac4$mean.Cor.over.CVs<=OneSEOfMax] <- NA_real_

gplots::heatmap.2(m.thresholded[-2,][c(3,1,4,2,5:9),], dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none') #used to fix unusual input value used to run the code

sum(alpha4$obj.ALPHA != 0)
# 436
sum(beta4$obj.BETA != 0)
# 248

#The highest lambdas with mean cor which is not more than 1 SE less than the max are:
