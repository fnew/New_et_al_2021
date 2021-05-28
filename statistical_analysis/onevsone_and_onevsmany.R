library(statmod)
library(glmnet)
#### simulate data... replace with real data when you can find it.. ####
n <- 10#0
p <- 100#0
q <- 200#0

X <- matrix(rnorm(n*p), nrow=n, ncol=p) # host
Y <- exp(matrix(rnorm(n*q), nrow=n, ncol=q)) # microbiome

#### one vs one ####

onevone_cors <- matrix(NA, nrow=ncol(X)/2, ncol=ncol(Y))
onevone_ps <- matrix(NA, nrow=ncol(X)/2, ncol=ncol(Y))
for (i in 1:nrow(onevone_cors)) {
  for (j in 1:ncol(onevone_cors)) {
    reg <- lm(Y[,j]~X[,2*i-1]+X[,2*i])#, family=tweedie(var.power=1.73, link.power=0))
    onevone_cors[i,j] <- with(summary(reg), r.squared)
    onevone_ps[i,j] <- with(summary(reg), r.squared)
  }
}
#onevone_cors[is.na(onevone_cors)] <- -1

orderstat <- Rfast::nth(onevone_cors, 100, descending = T)
entries_one <- which(onevone_cors >= orderstat, arr.ind = TRUE)
entries_one

#### one vs many ####

onevmany_cors <- matrix(NA, nrow=ncol(X)/2, ncol=ncol(Y))
for (j in 1:ncol(onevmany_cors)) {
  reg <- glmnet(X, Y[,j], nlambda=100, alpha=0)#, family=tweedie(var.power=1.73, link.power=0),) # use lasso here
  coefs <- coef(reg)[,ncol(coef(reg))][-1]
  for (k in 1:(length(coefs)/2)) {
    onevmany_cors[k,j] <- coefs[2*k-1]^2+coefs[2*k]^2
  }
  onevmany_cors[,j] <- sqrt(onevmany_cors[,j]/sum(onevmany_cors[,j]))
}
#onevmany_cors <- abs(onevmany_cors)
#onevmany_cors[is.na(onevmany_cors)] <- -1

orderstat <- Rfast::nth(onevmany_cors, 100, descending = T)
entries_many <- which(onevmany_cors >= orderstat, arr.ind = TRUE)
entries_many
