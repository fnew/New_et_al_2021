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