##########################################################################
##	Functions to perform SCCA 
##	Date: 2019
##
##	runCCA
##	get_cross_validated_penalty_parameters
##	get_enet
##	get_glasso
##	get_non_parallel_cv
##	get_parallel_cv
##	get_split_sets
##
##########################################################################


## Function to perform CCA ##
runCCA <- function(predictor,
                   predicted,
                   lam.eps = .01,
                   nlambda = 16,
                   penalization_x = "enet",
                   ridge_penalty_x = .95,
                   group_x = NA,
                   ALPHA = c(),
                   #nonzero_x = 1,
                   penalization_y = "enet",
                   ridge_penalty_y = .95,
                   group_y = NA,
                   BETA = c(),
                   #nonzero_y = 1,
                   #stop criteruim in a form of maximum iterations
                   max_iterations = 100, 
                   #stop criterium regarding to CRT
                   tolerance = 1*10^-20,
                   cross_validate = TRUE,
                   parallel_CV = TRUE,
                   nr_subsets = 10,
                   multiple_LV = FALSE,
                   nr_LVs = 1,
                   nr_cores = 1
){
  
  #print("SCCA started")
  ####################################################
  Y.mat <- as.matrix(predicted)
  #scale (SD = 1) and center (mean = 0) Y
  Yc <- scale(Y.mat)
  Yc <- Yc[,!colSums(!is.finite(Yc))]
  
  #print("Yc created")

  X.mat <- as.matrix(predictor)
  #scale (SD = 1) and center (mean = 0) X
  Xc <- scale(X.mat)
  Xc <- Xc[,!colSums(!is.finite(Xc))]
  
  #print("Xc created")
  
  n <- nrow(Xc)
  p <- ncol(Xc)
  q <- ncol(Yc)
  
  ########################################################
  Sigma12 <- crossprod(Xc, (Yc/n))
  #print("Sigma12 done")
  #Sigma1 <- crossprod(Xc, (Xc/n))
  #print("Sigma1 done")
  #Sigma2 <- crossprod(Yc, (Yc/n))
  #print("Sigma2 created")
  ########################################################
  if (length(ALPHA) == 0) {
    ALPHA <- rep(1,p)
    #for (j in 1:20) {
    #  APLHA <- t(Xc) %*% Xc %*% ALPHA
    #  ALPHA <- ALPHA / sqrt(sum(ALPHA^2))
    #}
  }
  ALPHA <- ALPHA/sqrt(as.numeric(crossprod(Xc %*% ALPHA/n, Xc %*% ALPHA)))
  #print("ALPHA made")
  if (length(BETA) == 0) {
    BETA <- rep(1,q)
    #for (j in 1:20) {
    #  BETA <- t(Yc) %*% Yc %*% BETA
    #  BETA <- BETA / sqrt(sum(BETA^2))
    #}
  }
  BETA <- BETA/sqrt(as.numeric(crossprod(Yc %*% BETA/n, Yc %*% BETA)))
  #print("BETA made")

  ###########################################################
  XtY <- Sigma12 %*% BETA
  YtX <- crossprod(Sigma12, ALPHA)
  
  XtYi <- max(abs(XtY))
  XtY2 <- sqrt(sum(XtY^2))
  YtXi <- max(abs(XtY))
  YtX2 <- sqrt(sum(YtX^2))
  #print("YtX2 made")


  lasso_penalty_x <- 1
  lasso_penalty_y <- 1
  grp_penalty_x <- 1
  grp_penalty_y <- 1
  
  if (penalization_x == "enet") {
    lasso_penalty_x_m <- (2/n)*(XtYi/ridge_penalty_x)
    lasso_penalty_x <- exp(seq(log(lasso_penalty_x_m), log(lam.eps*lasso_penalty_x_m), length.out = nlambda))
  }
  
  if (penalization_x == "glasso") {
    vec <- vector(length = length(XtY)/2)
    for (i in 1:(length(XtY)/2)) {
      my_seq <- c(2*i-1, 2*i)
      vec[i] <- sqrt(sum((XtY[my_seq])^2))
    }
    
    grp_penalty_x_m <- (sqrt(2)/n)*max(vec)
    grp_penalty_x <- exp(seq(log(grp_penalty_x_m), log(lam.eps*grp_penalty_x_m), length.out = nlambda))
  }
  
  if (penalization_y == "enet") {
    lasso_penalty_y_m <- (2/n)*(YtXi/ridge_penalty_y)
    lasso_penalty_y <- exp(seq(log(lasso_penalty_y_m), log(lam.eps*lasso_penalty_y_m), length.out = nlambda))
  }
  
  if (penalization_y == "glasso") {
    vec <- vector(length = length(XtY)/2)
    for (i in 1:(length(YtX)/2)) {
      my_seq <- c(2*i-1, 2*i)
      vec[i] <- sqrt(sum((YtX[my_seq])^2))
    }
    
    grp_penalty_y_m <- (sqrt(2)/n)*max(vec)    
    grp_penalty_y <- exp(seq(log(grp_penalty_y_m), log(lam.eps*grp_penalty_y_m), length.out = nlambda))
  }
  
  obj <- sCCA(predictor = predictor,
              predicted = predicted,
              penalization_x = penalization_x,
              lasso_penalty_x = lasso_penalty_x,
              grp_penalty_x = grp_penalty_x,
              ridge_penalty_x = ridge_penalty_x,
              group_x = group_x,
              ALPHA = ALPHA,
              #nonzero_x = 1,
              penalization_y = penalization_y,
              lasso_penalty_y = lasso_penalty_y,
              grp_penalty_y = grp_penalty_y,
              ridge_penalty_y = ridge_penalty_y,
              group_y = group_y,
              BETA = BETA,
              #nonzero_y = 1,
              #stop criteruim in a form of maximum iterations
              max_iterations = max_iterations,
              #stop criterium regarding to CRT
              tolerance = tolerance,
              cross_validate = cross_validate,
              parallel_CV = parallel_CV,
              nr_subsets = nr_subsets,
              multiple_LV = multiple_LV,
              nr_LVs = nr_LVs,
              nr_cores = nr_cores)
  
  return(obj)
}





### FUNCTION TO CROSS VALIDATE PENALTY PARAMETERS ###
get_cross_validated_penalty_parameters <- function(predictor,
                                                   predicted,
                                                   penalization_x, #
                                                   lasso_penalty_x,
                                                   grp_penalty_x,
                                                   ridge_penalty_x,
                                                   group_x,
                                                   #nonzero_x,
                                                   penalization_y, #
                                                   lasso_penalty_y,
                                                   grp_penalty_y,
                                                   ridge_penalty_y,
                                                   group_y,
                                                   #nonzero_y,
                                                   nr_subsets, #
                                                   max_iterations,
                                                   tolerance,
                                                   parallel_CV = parallel_CV,
                                                   nr_cores #added
                                                   ){
  #print("Started get_cross_validated_penalty_parameters ")
  shuffled <-  get_split_sets(X = predictor,
                              Y = predicted,
                              nr_subsets = nr_subsets)

  X.sampled     <-   shuffled$X.sampled
  Y.sampled     <-   shuffled$Y.sampled
  label         <-   shuffled$labels

  if (parallel_CV){
    #print("doing parallel") #
    cv_results <- get_parallel_cv(X = X.sampled,
                                  Y = Y.sampled,
                                  penalization_x = penalization_x, #####
                                  lasso_penalty_x = lasso_penalty_x,
                                  grp_penalty_x = grp_penalty_x,
                                  ridge_penalty_x = ridge_penalty_x,
                                  group_x = group_x,
                                  #nonzero_x = nonzero_x,
                                  penalization_y = penalization_y, #####
                                  lasso_penalty_y = lasso_penalty_y,
                                  grp_penalty_y = grp_penalty_y,
                                  ridge_penalty_y = ridge_penalty_y,
                                  group_y = group_y,
                                  #nonzero_y = nonzero_y,
                                  label = label, #######################
                                  #penalization = penalization,
                                  max_iterations = max_iterations,
                                  tolerance = tolerance,
                                  parallel_CV = parallel_CV, ##BB added
                                  nr_cores) ##BB added

  } else {
    #print("not doing parallel") #
    cv_results <- get_non_parallel_cv(X = X.sampled,
                                      Y = Y.sampled,
                                      penalization_x = penalization_x, #####
                                      lasso_penalty_x = lasso_penalty_x,
                                      grp_penalty_x = grp_penalty_x,
                                      ridge_penalty_x = ridge_penalty_x,
                                      group_x = group_x,
                                      #nonzero_x = nonzero_x,
                                      penalization_y = penalization_y, #####
                                      lasso_penalty_y = lasso_penalty_y,
                                      grp_penalty_y = grp_penalty_y,
                                      ridge_penalty_y = ridge_penalty_y,
                                      group_y = group_y,
                                      #nonzero_y = nonzero_y,
                                      label = label,
                                      #penalization = penalization,
                                      max_iterations = max_iterations,
                                      tolerance = tolerance,
                                      parallel_CV = parallel_CV) ##BB added

  }
  
  #print("finished cross validation") #  
  
  a = cv_results$mean_abs_cors[,7] #BB: changed from 3 to 7
  #print(a) #

  best_values <- cv_results$mean_abs_cors[which.max(a),]
  #print(best_values)
  
  best_ridge_y <- best_values[1]
  best_lasso_y <- best_values[2]
  best_group_y <- best_values[3]
  best_ridge_x <- best_values[4]
  best_lasso_x <- best_values[5]
  best_group_x <- best_values[6]
  #best_ridge <- best_values[1]
  #best_nonzero <- best_values[2]


  #**********************
  result <- list(
                 abs_cors = cv_results$abs_cors,
                 mean_abs_cors = cv_results$mean_abs_cors,
                 stime = cv_results$stime,
                 iterations_m = cv_results$iterations_m,
                 #best_ridge = best_ridge,
                 #best_nonzero = best_nonzero
                 best_ridge_x = best_ridge_x,
                 best_group_x = best_group_x,
                 best_lasso_x = best_lasso_x,
                 best_ridge_y = best_ridge_y,
                 best_group_y = best_group_y,
                 best_lasso_y = best_lasso_y
  )
  #class(result) = "CVsRDA"
  result
}


#'@S3method print CVsRDA
print.CVsRDA <- function(x, ...)
{
  cat("Cross validation (CV) for sRDA or sCCA", "\n")
  cat(rep("-", 5), sep="")
  cat("\n   NAME                ", "DESCRIPTION")
  cat("\n1  $abs_cors           ", "sum of absolute correlations per k-th fold CV")
  cat("\n2  $mean_abs_cors      ", "mean absolute correlation per CV")
  cat("\n3  $best_ridge         ", "best ridge parameter selected for the model")
  cat("\n4  $best_nonzero       ", "best number of nonzero alpha weights selected")
  cat("\n5  $stime              ", "time elapsed in seconds")
  cat("\n6  $iterations_m       ", "number of iterations ran before convergence")
  cat("\n")
  cat(rep("-", 45), sep="")
  invisible(x)
}


### FUNCTION FOR ELASTIC NET PENALIZATION FOR MLR ###
  #input:     X         (matrix - independent variables),
  #           y         (matrix - dependent variables),
  #           lambda    (integer - factor for Ridge penalty)
  #           nonzero   (integer - factor for LASSO penalty)
  
  #output:    tmp       vector    peanalized weights after regressing y on x
get_enet = function(X, y, lambda, alpha){
 	as.numeric(glmnet::glmnet(X, y, lambda = lambda, alpha = alpha)$beta)
}


### FUNCTION FOR GROUP LASSO PENALIZATION ###
get_glasso <- function(X, y, group, lambda) {
  beta <- c(gglasso::gglasso(X, y, group = group, lambda = lambda)$beta)
  names(beta) <- colnames(X)
  return(beta)
}



### CROSS VALIDATION FOR FINDING OPTIMAL PENALIZATION PARAMS ###
get_non_parallel_cv <- function(X,
                                Y, ####################
                                penalization_x, #######
                                lasso_penalty_x,
                                grp_penalty_x,
                                ridge_penalty_x,
                                group_x,
                                #nonzero_x,
                                penalization_y, #######
                                lasso_penalty_y,
                                grp_penalty_y,
                                ridge_penalty_y,
                                group_y,
                                #nonzero_y,
                                label, ################
                                #penalization,
                                max_iterations,
                                tolerance,
                                parallel_CV){ #BB added
  nr_subsets      <-    length(unique(label))
  abs_cors        <-    c()
  iterations_m    <-    c()
  kth_fold        <-    0
  
  #print("got into get_non_parallel_cv") #
  ii <- 0 
  
  #Measure time
  stime <- system.time({
    for (l2x in 1:length(grp_penalty_x)){
    for (l1x in 1:length(lasso_penalty_x)){
    for (l3x in 1:length(ridge_penalty_x)){
      for (l2y in 1:length(grp_penalty_y)){
      for (l1y in 1:length(lasso_penalty_y)){
      for (l3y in 1:length(ridge_penalty_y)){
        
        ii <- ii+1 
	      print(paste("##########################Doing tuning param combo #", ii, '')) 

        sub_abs_cor       <- c() ##BB: moved here from 2 lines up
        sub_iterations_m  <- c() ##BB: moved here from 2 lines up
        
        ALPHA.old <- c() 
	      BETA.old <- c()

        for (i in 1:nr_subsets){
          print(paste("doing subset number", i))

          kth_fold <- kth_fold + 1 ##BB: moved here from 2 lines up
          
          X.train   <- X[label!=i,]
          X.test    <- X[label==i,]
          
          Y.train   <- Y[label!=i,]
          Y.test    <- Y[label==i,]
          
          sub_results <- sCCA(predictor = X.train, #removed [[i]]
                                   predicted = Y.train,
                                   penalization_x = penalization_x, #######
                                   lasso_penalty_x = lasso_penalty_x[l1x],
                                   grp_penalty_x = grp_penalty_x[l2x],
                                   ridge_penalty_x = ridge_penalty_x[l3x],
                                   group_x = group_x,
                                   ALPHA = ALPHA.old,
                                   #nonzero_x,
                                   penalization_y = penalization_y, #######
                                   lasso_penalty_y = lasso_penalty_y[l1y],
                                   grp_penalty_y = grp_penalty_y[l2y],
                                   ridge_penalty_y = ridge_penalty_y[l3y],
                                   group_y = group_y,
                                   BETA = BETA.old,
                                   #nonzero_y,
                                   #####penalization = penalization, #########
                                   max_iterations = max_iterations,
                                   tolerance = tolerance,
                                   cross_validate = FALSE,
                                   parallel_CV = parallel_CV) #BB: added
          
          ALPHA.old <- sub_results$ALPHA # removed [[i]] 
	        BETA.old <- sub_results$BETA # removed [[i]]

	        
	        #print(paste0('dim of X.test:', dim(X.test)))
	        #print(paste0('length of ALPHA.old: ', length(ALPHA.old)))
	        
	        
	        
          XI.test <- scale(X.test) %*% ALPHA.old
          
          #Divide with dim(Y.train)[2], exclude NA's from sum
          sub_abs_cor[[i]]            <- sum(abs(cor(XI.test,Y.test)),na.rm = T) #atomic vector
          sub_iterations_m[[i]]       <- sub_results$nr_iterations

        } #End of subset for loop

        abs_cors        <- cbind(abs_cors, sub_abs_cor) #matrix
        iterations_m    <- cbind(iterations_m, sub_iterations_m)

      } #End of tuning params loop
      }
      }
    }
    }
    } #End of tuning params loop
  })[3] #End of measure time
  
  #print("Exited tuning looping")

  labels_tuning <- expand.grid(ridge_penalty_y, lasso_penalty_y, grp_penalty_y, 
                               ridge_penalty_x, lasso_penalty_x, grp_penalty_x) ##BB: added

  mean_abs_cors <- c()
  for (i in 1:nrow(labels_tuning)){
    mean_abs_cors <- rbind(mean_abs_cors, #BB: this is the existing rows
                           c(labels_tuning[i,], mean(abs_cors[,i])))
  }

  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Y Ridge Penalty", "Y Lasso Penalty", "Y Group Penalty",
                                 "X Ridge Penalty", "X Lasso Penalty", "X Group Penalty",
                                 "mean Cor over CVs")

  ## Return section ##
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime,
                           iterations_m = iterations_m)

  #print("Completed get_non_parallel")
  result
}


### PARALLEL N-FOLD CROSS VALIDATION ###
get_parallel_cv <- function(X,
                            Y, ####################
                            penalization_x, ####### BB
                            lasso_penalty_x,
                            grp_penalty_x,
                            ridge_penalty_x,
                            group_x,
                            #nonzero_x,
                            penalization_y, #######
                            lasso_penalty_y,
                            grp_penalty_y,
                            ridge_penalty_y,
                            group_y,
                            #nonzero_y, BB
                            label, ################
                            #penalization, BB
                            max_iterations,
                            tolerance,
                            parallel_CV,
                            nr_cores){

  nr_subsets    <-    length(unique(label))
  #print("Doing parallel CV")

  ## PARALLELIZATION ##
  cl    <-    parallel::makeForkCluster(nr_cores) ###BB: can you even use all your cores at once??
  doParallel::registerDoParallel(cl)
  
  #Measure time
  stime <- system.time({
    #Get the alphas with parallel computing
    x <- foreach(grp_pen_x = grp_penalty_x, .export= c("sCCA", "get_glasso", "get_enet"), .combine = c) %dopar% { 
      abs_cors          <- c()
      iterations_m      <- c()
      
      #for (l2x in 1:length(grp_penalty_x)){
      for (l1x in 1:length(lasso_penalty_x)){
      for (l3x in 1:length(ridge_penalty_x)){
      for (l2y in 1:length(grp_penalty_y)){
      for (l1y in 1:length(lasso_penalty_y)){
      for (l3y in 1:length(ridge_penalty_y)){
        
        sub_abs_cor       <- c()
	      sub_iterations_m  <- c() 
        
        ALPHA.old <- c() 
	      BETA.old <- c()

        for (i in 1:nr_subsets){
          X.train   <- X[label!=i,]
          X.test    <- X[label==i,]

          Y.train   <- Y[label!=i,]
          Y.test    <- Y[label==i,]

          sub_results <- sCCA(predictor = X.train,
                              predicted = Y.train,
                              penalization_x = penalization_x, #######
                              lasso_penalty_x = lasso_penalty_x[l1x],
                              grp_penalty_x = grp_pen_x,
                              ridge_penalty_x = ridge_penalty_x[l3x],
                              group_x = group_x,
                              ALPHA = ALPHA.old,
                              #nonzero_x,
                              penalization_y = penalization_y, #######
                              lasso_penalty_y = lasso_penalty_y[l1y],
                              grp_penalty_y = grp_penalty_y[l2y],
                              ridge_penalty_y = ridge_penalty_y[l3y],
                              group_y = group_y,
                              BETA = BETA.old,
                              #nonzero_y,
                              #####penalization = penalization, #########
                              max_iterations = max_iterations,
                              tolerance = tolerance,
                              cross_validate = FALSE)
                              #parallel_CV = FALSE) #BB: added #Copypaste from nonparallel #not sure why this is added here
          
          ALPHA.old <- sub_results$ALPHA 
	        BETA.old <- sub_results$BETA 

          XI.test = scale(X.test) %*% ALPHA.old

          sub_iterations_m[[i]] <- sub_results$nr_iterations ##BB: NB this is an atomic vector (ult. length = nr_subsets)
          sub_abs_cor[[i]] <- sum(abs(cor(XI.test,Y.test)), na.rm = TRUE) ##BB: why sum here?
        } #End of subset for loop
        
        abs_cors        <- cbind(abs_cors, sub_abs_cor) ##BB: NB this is a matrix (ult. nr_subsets x (#tuningparamsbesidesgrpx))
        iterations_m    <- cbind(iterations_m, sub_iterations_m)
      } #End of all other params loop
      }
      }
      }
      }
      list(abs_cors, iterations_m)
    } #End of lambda for loop
  })[3] #End of measure time

  #Stop cluster
  stopCluster(cl) 
  
  #mine out the iterations from vector x 
  abs_cors <- c()
  iterations_m <- c()
  j <- 0
  for (i in 1:length(grp_penalty_x)) {
    j <- j+1
    abs_cors <- cbind(abs_cors, x[[j]])
    j <- j+1
    iterations_m <- cbind(iterations_m, x[[j]])
  }
  
  #Arrange like in get_non_parallel 
  labels_tuning <- expand.grid(ridge_penalty_y, lasso_penalty_y, grp_penalty_y, 
                               ridge_penalty_x, lasso_penalty_x, grp_penalty_x) ##BB: added

  mean_abs_cors <- c()
  for (i in 1:nrow(labels_tuning)){
    mean_abs_cors <- rbind(mean_abs_cors, #BB: this is the existing rows
                           c(labels_tuning[i,], mean(abs_cors[,i])))
  }

  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Y Ridge Penalty", "Y Lasso Penalty", "Y Group Penalty",
                                 "X Ridge Penalty", "X Lasso Penalty", "X Group Penalty",
                                 "mean Cor over CVs")
  
  #Return section
  result        <-    list(abs_cors = abs_cors,
                           mean_abs_cors = mean_abs_cors,
                           stime = stime)
  result
}


### FUNCTION TO SPLIT UP DATA MATRIX INTO n SUBSETS ###
get_split_sets <- function(X, Y, nr_subsets){
#
#input:     X             n*p matrix      - data matrix to be split
#           Y             n*q matrix      - data matrix to be split
#           nr_subsets    integer         - nr of subsets the data
#                                           matrix will be split to
#
#
#output:    label         vector with n_subset unique elements with length of n
#
#           X.sampled,    n*p matrix      shuffled dataset of the origninal data
#           Y.sampled,    n*q matrix      shuffled dataset of the origninal data

#print("Running get_split_sets")

  n <- dim(X)[1]

  #re-sample data rows
  splitting_dimensions <- sample(1:n,n)

  X.sampled <- X[splitting_dimensions,]
  Y.sampled <- Y[splitting_dimensions,]

  #calculate how many rows are left over
  leftover = n %% nr_subsets

  rep_eat = (n-leftover)/nr_subsets

  #repeat sequence nr_subsets
  labels = rep(1:nr_subsets , each=rep_eat)

  if(leftover!=0)
  {
    labels = c(labels, (1:nr_subsets)[1:leftover])
  }

  #Return section
  result <- list(X.sampled = X.sampled,
                 Y.sampled = Y.sampled,
                 labels = labels)
  result
}


sCCA <- function(predictor,
                 predicted,
                 penalization_x = "enet",
                 lasso_penalty_x = 1,
                 grp_penalty_x = 1,
                 ridge_penalty_x = 1,
                 group_x = NA,
                 ALPHA = c(),
                 #nonzero_x = 1,
                 penalization_y = "enet",
                 lasso_penalty_y = 1,
                 grp_penalty_y = 1,
                 ridge_penalty_y = 1,
                 group_y = NA,
                 BETA = c(),
                 #nonzero_y = 1,
                 #stop criteruim in a form of maximum iterations
                 max_iterations = 100,
                 #stop criterium regarding to CRT
                 tolerance = 1*10^-20,
                 cross_validate = FALSE,
                 parallel_CV = TRUE,
                 nr_subsets = 10,
                 multiple_LV = FALSE,
                 nr_LVs = 1,
                 nr_cores = 1
){
  #print("Starting SCCA")
  #****************************************************************************#
  #0. ancillary functions
  #*********************
  #check for dependencies
  if (!requireNamespace("elasticnet", quietly = TRUE)) {
    stop("The elasticnet package is needed for this function to work.
         Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop("The stats package is needed for this function to work.
         Please install it.",
         call. = FALSE)
  }

  #if ust penalization, change ridge penalties to Inf
  if (penalization_x == "ust"){
    ridge_penalty_x <- rep(Inf, length(ridge_penalty))
  }
  if (penalization_y == "ust"){
    ridge_penalty_y <- rep(Inf, length(ridge_penalty))
  }

  #check for dependencies for parallel computing
  if (parallel_CV == TRUE){

    #check for dependencies
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The parallel package is needed for this function to work.
           Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("The foreach package is needed for this function to work.
           Please install it.",
           call. = FALSE)
    }

    #check for cores
    nr_cores      <-   nr_cores
    if(nr_cores){
      cat("Cross valdiation function on multiple cores is called,
          number of cores to be used on this work station:", nr_cores, "\n")
    }else{
      stop("R wasn't able to identify multiple cores on this work station.
           Please make sure your work station is compatible with
           parallel computing.",
           call. = FALSE)
    }

  }

  #Check input values for the function
  if( dim(predictor)[1] != dim(predicted)[1]){
    stop("The sample size for the predictor and predicted data sets is
         different. Please use data sets with equal sample size.",
         call. = FALSE)
  }

  #Check for cross validation *************************************************#

#  if (multiple_LV){
#
#    result <- get_multiple_LVs(X = predictor,
#                               Y = predicted,
#                               penalization = penalization,
#                               lambda = ridge_penalty,
#                               nonzero = nonzero,
#                               nr_latent = nr_LVs,
#                               stop_criterium = tolerance,
#                               max_iterations = max_iterations,
#                               cross_validate = cross_validate)
#
#    return(result)
#
#
#  }
  
  CV_results <- "Cross validation function is not called"

  if (cross_validate){

    if (penalization_x=="none" && penalization_y=="none") {
      stop("You do not need cross validation for a non penalized model
           (pelase change the penalization or cross_validate variable)",
           call. = FALSE)
    }
    #print("entering get_cross_validated_penalty_parameters") ########################################
    CV_results <-
      get_cross_validated_penalty_parameters(predictor = predictor,
                                             predicted = predicted,
                                             penalization_x = penalization_x,
                                             lasso_penalty_x = lasso_penalty_x,
                                             grp_penalty_x = grp_penalty_x,
                                             ridge_penalty_x = ridge_penalty_x,
                                             group_x = group_x,
                                             penalization_y = penalization_y,
                                             lasso_penalty_y = lasso_penalty_y,
                                             grp_penalty_y = grp_penalty_y,
                                             ridge_penalty_y = ridge_penalty_y,
                                             group_y = group_y,
                                             nr_subsets = nr_subsets,
                                             max_iterations = max_iterations,
                                             tolerance = tolerance,
                                             parallel_CV = parallel_CV,
                                             nr_cores = nr_cores
      )
    
    ridge_penalty_x <- CV_results$best_ridge_x[[1]] #BB: having to unlist
    grp_penalty_x <- CV_results$best_group_x[[1]]
    lasso_penalty_x <- CV_results$best_lasso_x[[1]]
    ridge_penalty_y <- CV_results$best_ridge_y[[1]]
    group_penalty_y <- CV_results$best_group_y[[1]]
    lasso_penalty_y <- CV_results$best_lasso_y[[1]]
    

  } else {
    # if no cross validation and more then 1 ridge/non zero penalty, give error
    if (length(ridge_penalty_x)>1 || length(grp_penalty_x)>1 || length(lasso_penalty_x) > 1 || 
        length(ridge_penalty_y)>1 || length(grp_penalty_y)>1 || length(lasso_penalty_y) > 1) { #BB:edited to include all
      stop("Multiple tuning parameters where only one is
           needed (pelase set cross_validate = TRUE if you would like to run
           cross validation)",
           call. = FALSE)
    }

  }

  #Algorithm for peanalization*************************************************#

  #****************************************************************************#
  #1. Preperation of the data
  #*************************#

  Y.mat <- as.matrix(predicted)
  #scale (SD = 1) and center (mean = 0) Y
  Yc <- scale(Y.mat)
  Yc <- Yc[,!colSums(!is.finite(Yc))]

  X.mat <- as.matrix(predictor)
  #scale (SD = 1) and center (mean = 0) X
  Xcr <- scale(X.mat)
  Xcr <- Xcr[,!colSums(!is.finite(Xcr))]

  #variable length of X and Y
  p <- ncol(Xcr)
  q <- ncol(Yc)

  #   For initalization, let:
  #
  #         ETA^(0) = y_1 + y_2 + .... + y_q = Y*BETA_hat^(0) ( = Y_i)
  #           were BETA^(0) = [1,1....,1]
  # an y-weight column vector of q elements
  if (length(BETA)==0) {
    BETA = matrix( rep(1, q), nrow = q, ncol = 1, byrow = TRUE)
   #print("entered, ie had length(beta) = 0") ########################################################
  }
  # an x-weight column vector of p elements
  if (length(ALPHA)==0) {
    ALPHA = matrix( rep(1, p), nrow = p, ncol = 1, byrow = TRUE)
  }

  #if(cross_validate) {
  #  print("here in main function, after tuning") ########################################################
  #}
  #print(length(ALPHA)) ########################################################
  #print(length(BETA)) ########################################################
  
  #****************************************************************************#
  #2. Iterative loop until convergence
  #*************************#
  CRTs          <- c()
  sum_abs_Betas <- c()
  Nr_iterations = 0         #iteration counter

  WeContinnue = TRUE        #T/F value for stop criterium based on CRT values
  #in last two iteration
  CRT = 1                   #convergance measure between alpha and beta

  #print("started iterating in sCCA") ###############################################################
  print(colnames(Yc)[1:10])
  while(CRT > tolerance && WeContinnue && Nr_iterations < max_iterations) {

    ETA = Yc %*% BETA
    ETA = scale(ETA)

    XI = Xcr %*% ALPHA
    XI = scale(XI)

    #**************************************************************************#
    #Chose between generic RDA, SRDA with ENET or SRDA wiht UST****************#

    # For the value of ahat^(0) (ALPHA) and hence XIhat^(0),
    # regress ETAhat^(0) on X jointly to get ALPHA,
    # ALPH_0 = (X'X)^-1 (X' ETA)

    ALPH_0 <- switch(penalization_x,
                     #lm without penalization
                     #ALPH_0 = as.matrix(lm(ETA ~ 0+Xcr)$coefficients)
                     "none" = {
                       solve(t(Xcr) %*% Xcr) %*% t(Xcr) %*% ETA
                     },
                     
                     #CALCULATE WITH ELASTIC NET
                     #ALPH_0 = as.matrix(calculateVectorEnet(Xcr,ETA,1,30))
                     "enet" = {
                       as.matrix(get_enet(Xcr,ETA,lasso_penalty_x, ridge_penalty_x))
                     },
                     
                     #CALCULATE WITH GROUP LASSO
                     "glasso" = {
                       as.matrix(get_glasso(Xcr, ETA, group_x, grp_penalty_x))
                     },

                     #CALCULATE WITH UST
                     #ALPH_0 = as.matrix(calculateVectorEnet(Xcr,ETA,1,30))
##                     "ust" = {
##                       as.matrix(get_ust(Xcr,ETA,nonzero))
##                     },
                     {
                       stop("Please choose a valid penalization method
                            (i.e. 'ust', 'enet' or 'none')",
                            call. = FALSE)

                     }

                       )

    #***************************************************************************
    #         compute the value of XI^(1):
    #         XIhat^(1) = SUM_running to p where t=1 ( ahat_t^(0) *x_t )

    XI = Xcr %*% ALPH_0

    #         and normalize XIhat^(1) such that
    #           t(XIhat^(1))*XIhat^(1) = t(ahat^(0)) t(X)*X*ahat^(0) = 1
    #           t(XI)%*%XI = t(ALPH_0) %*% t(Xcr) %*% Xcr %*% ALPH_0 = 1
    #           that is its variance is 1
    XI = scale(XI)
    
    if (is.nan(XI[1])) {
      #print("XI IS NAN")
      XI <- rep(0, length(XI))
      ETA <- rep(0, length(ETA))
      ALPHA <- matrix(rep(0, length(ALPHA)), ncol=1)
      BETA <- matrix(rep(0, length(BETA)), ncol=1)
      break
    }

    #         For the value BETAhat^(1) and hence ETAhat^(1),
    #           regress y1,y2 ... yq separately on XIhat^(1),
    #
    #           y1 = BETA_1 * XIhat^(1) + Epsilon_1
    #           .
    #           yq = BETA_q * XIhat^(1) + Epsilon_q

    BETA_0 <- switch(penalization_y,
                     #lm without penalization
                     "none" = {
                       solve(t(Yc) %*% Yc) %*% t(Yc) %*% XI
                     },

                     #CALCULATE WITH ELASTIC NET
                     "enet" = {
                       as.matrix(get_enet(Yc,XI,lasso_penalty_y,ridge_penalty_y))
                     },
                     
                     #CALCULATE WITH GROUP LASSO
                     "glasso" = {
                       as.matrix(get_glasso(Yc, XI, group_y, grp_penalty_y))
                     },

##                     #CALCULATE WITH UST
##                     "ust" = {
##                       as.matrix(get_ust(Yc,XI,nonzero))
##                     },
                     {
                       stop("Please choose a valid penalization method
                            (i.e. 'ust', 'enet' or 'none')",
                            call. = FALSE)

                     }

                       )

    #BETA
    #         compute the vaule of ETAhat^(1),
    #
    #           ETAhat^(1) = SUM_running to q where k=1 ( BETAhat_k^(1) * y_k)
    ETA = Yc %*% BETA_0

    ETA = scale(ETA)
    
    if (is.nan(ETA[1])) {
      #print("ETA IS NAN")
      XI <- rep(0, length(XI))
      ETA <- rep(0, length(ETA))
      ALPHA <- matrix(rep(0, length(ALPHA)), ncol=1)
      BETA <- matrix(rep(0, length(BETA)), ncol=1)
      break
    }

    #check if ETA needed to be scaled here <- nope
    #ETA = scale(ETA)


    #Calculate convergence of Alphas and Betas
    CRT = sum((ALPHA - ALPH_0)^2, (BETA - BETA_0)^2)/(length(ALPHA) + length(BETA)); #BB: added division by lengths


    ALPHA=            ALPH_0
    BETA =            BETA_0

    #Check if last two iterations CR converges*********************************#
    Nr_iterations                   =   Nr_iterations + 1
    CRTs[[Nr_iterations]]           =   CRT
    sum_abs_Betas[[Nr_iterations]]  =   sum(abs(BETA))

    if (Nr_iterations>1){

      stop_condition <- abs(CRTs[[Nr_iterations]] - CRTs[[Nr_iterations-1]])
      stop_criterium <- 1 * 10^-6

      if (stop_condition < stop_criterium){
        WeContinnue <- FALSE
      }

    }#END Check if last two iterations CR converges*************************#

  }# End of main loop
  
  print(Nr_iterations)

  #REDN = REDUNDANCY INDEX: redundancy of y variables
  #  REDN = (sum(BETA^2))/4
  Red_index = (sum(BETA^2))/4

  #************************

  #  t(BETA) %*% t(XYMA) %*% solve(XXMA) %*% XYMA %*% BETA

  #  PATH


  #***************************

  #use this later for second latent variable
  #SOLVE_XIXI = solve(t(XI)%*%XI) %*% t(XI)       #BB: REMOVED

  rownames(ALPHA) <- colnames(Xcr)
  rownames(BETA) <- colnames(Yc)

  result <- list(XI = XI,
                 ETA = ETA,
                 ALPHA = ALPHA,
                 BETA= BETA,
                 nr_iterations = Nr_iterations,
                 #inverse_of_XIXI = SOLVE_XIXI,  #BB: REMOVED
                 iterations_crts = CRTs,
                 sum_absolute_betas = sum_abs_Betas, #############
                 penalization_x = penalization_x,
                 lasso_penalty_x = lasso_penalty_x,
                 grp_penalty_x = grp_penalty_x,
                 ridge_penalty_x = ridge_penalty_x,
                 penalization_y = penalization_y,
                 lasso_penalty_y = lasso_penalty_y,
                 grp_penalty_y = grp_penalty_y,
                 ridge_penalty_y = ridge_penalty_y, #################
                 nr_latent_variables = nr_LVs,
                 CV_results = CV_results
  )

  class(result) = "sRDA"
  return(result)
}
