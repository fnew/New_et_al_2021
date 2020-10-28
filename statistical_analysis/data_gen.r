setwd('/Users/benbaer/Box Sync/school/research/felicia/cca_debugging/')
#### reading in functions ####
generate_data <- function(w1, w2, n) {
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  
  X <- c()
  Y <- c()
  for (i in 1:n) {
    X <- rbind(X, w1*z1[i] + sample(w1)*z2[i] + rnorm(length(w1)))
    Y <- rbind(Y, w2*z1[i] + sample(w2)*z2[i] + rnorm(length(w2)))
  }
  list(X, Y)
} # the signal strength could be the norms of w1 and w2, since when those are 0 there's no relationship between X and Y

nmlize <- function(a) a/sqrt(sum(a^2))

source('refactored_cca_functions.R')

scale_back <- function(ll) {
  alpha <- ll[[1]]
  beta <- ll[[2]]
  
  list('alpha'=nmlize(alpha/attributes(scale(X))$'scaled:scale'),
       'beta'=nmlize(beta/attributes(scale(Y))$'scaled:scale'))
}
sb_x <- function(alpha) {
  nmlize(alpha/attributes(scale(X))$'scaled:scale')
}
sb_y <- function(beta) {
  nmlize(beta/attributes(scale(Y))$'scaled:scale')
}
get_list <- function(out) list('alpha'=out$xcoef[,1], 'beta'=out$ycoef[,1])

#### testing that generate data recovers correctly on low dim data for cca ####
alpha <- nmlize(c(-10, -5, 0, 7)) #rpois(5,1)
beta <- nmlize(c(-10, 0, 15)) #rpois(7,3)
n <- 1000

gen_out <- generate_data(alpha, beta, n)
Y <- gen_out[[2]]
X <- gen_out[[1]]