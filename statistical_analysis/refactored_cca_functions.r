#### scca w tuning ####
scca_w_tuning <- function(Xdef, Ydef, n, p, q, 
                     lam_eps, nlambda,
                     penalization_x,relpenalty_x,group_x,alpha,
                     penalization_y,relpenalty_y,group_y,beta,
                     max_iterations,tolerance,
                     cross_validate,num_folds,parallel_CV,num_cores,comp=0) {
  ## initializing coefficients
  if (length(alpha) == 0) alpha <- rep(1,p)
  xi <- Xdef %*% alpha
  alpha <- alpha/sqrt(sum(xi^2)*n)
  
  if (length(beta) == 0) beta <- rep(1,q)
  eta <- Ydef %*% beta
  beta <- beta/sqrt(sum(eta^2)*n)
  
  # do cross validation if requested
  if (cross_validate) {
    ## get useful terms
    Sigma12 <- crossprod(Xdef, Ydef/n)
    
    XtY <- Sigma12 %*% beta
    YtX <- crossprod(Sigma12, alpha)
    
    overallpenalty_x <- 1
    overallpenalty_y <- 1
    
    ## define grid
    if (penalization_x == "enet") overallpenalty_x <- get_enet_grid(XtY, relpenalty_x, nlambda, lam_eps)
    if (penalization_y == "enet") overallpenalty_y <- get_enet_grid(YtX, relpenalty_y, nlambda, lam_eps)
    
    if (penalization_x == "glasso") overallpenalty_x <- get_glasso_grid(XtY, group_x, nlambda, lam_eps)
    if (penalization_y == "glasso") overallpenalty_y <- get_glasso_grid(YtX, group_y, nlambda, lam_eps)
    
    #print(overallpenalty_x)
    #print(overallpenalty_y)
    
    ## do cross validation
    CV_results <- get_cross_validated_penalty_parameters(
      X = Xdef, Y = Ydef,
      penalization_x = penalization_x, penalization_y = penalization_y,
      overallpenalty_x = overallpenalty_x, overallpenalty_y = overallpenalty_y,
      relpenalty_x = relpenalty_x, relpenalty_y = relpenalty_y,
      group_x = group_x, group_y = group_y,
      max_iterations = max_iterations, tolerance = tolerance,
      num_folds = num_folds, parallel_CV = parallel_CV, num_cores = num_cores
    )
    
    ## get output
    overallpenalty_x <- CV_results$best_overall_x
    overallpenalty_y <- CV_results$best_overall_y
  } else {
    overallpenalty_x <- overallpenalty_y <- .01 # arbitrarily chosen
    CV_results <- list(message="cross validation was not called")
  }
  
  # do algorithm
  out <- cca_algorithm(X=Xdef,Y=Ydef,
                       alpha=alpha,beta=beta,
                       penalization_x=penalization_x,penalization_y=penalization_y,
                       overallpenalty_x=overallpenalty_x, overallpenalty_y=overallpenalty_y, 
                       relpenalty_x=relpenalty_x, relpenalty_y=relpenalty_y, 
                       group_x=group_x, group_y=group_y, 
                       tolerance = tolerance, max_iterations = max_iterations
  )
  c(out, CV_results)
}

#### scca front end to do the data normalization, initing coefs, tuning parameter grid finding, and calling cross validation ####

#lam_eps=.01;nlambda=2;penalization_x="enet";relpenalty_x=.9;group_x=NA;alpha=c();penalization_y="enet";
#relpenalty_y=.9;group_y=NA;beta=c();max_iterations=100;tolerance=10^-20;cross_validate=TRUE;
#num_folds=5;parallel_CV=FALSE;num_cores=1;num_components=1;
scca_front <- function(X,Y,
                       lam_eps = .01, nlambda = 12,
                       penalization_x = "enet", relpenalty_x = .9, group_x = NA,
                       alpha = c(),
                       penalization_y = "enet", relpenalty_y = .9, group_y = NA,
                       beta = c(),
                       max_iterations = 100, tolerance = 10^-20,
                       cross_validate= TRUE, num_folds = 5, parallel_CV = FALSE, num_cores = 1,
                       num_components=1
) {
  # preprocessing data
  ## scaling data views and defining sizes
  Y <- as.matrix(Y)
  Yc <- scale(Y)
  Yc <- Yc[,!colSums(!is.finite(Yc))]
  
  X <- as.matrix(X)
  Xc <- scale(X)
  Xc <- Xc[,!colSums(!is.finite(Xc))]
  
  n <- nrow(Xc)
  p <- ncol(Xc)
  q <- ncol(Yc)
  
  # start looping
  list_comps <- list()
  A <- matrix(nrow=p,ncol=0)
  B <- matrix(nrow=q,ncol=0)
  U_comps <- matrix(nrow=n, ncol=0)
  V_comps <- matrix(nrow=n, ncol=0)
  
  Xdef <- Xc
  Ydef <- Yc
  for (comp in 1:num_components) {
    print(paste0('################Running component number ', comp))
    # get output from cca_alg after tunining and initing

    if (comp > 1) {
      alpha <- c()
      beta <- c()
    }
    
    out <- scca_w_tuning(Xdef, Ydef, n, p, q,
                         lam_eps=lam_eps, nlambda=nlambda,
                         penalization_x=penalization_x,relpenalty_x=relpenalty_x,group_x=group_x,alpha=alpha,
                         penalization_y=penalization_y,relpenalty_y=relpenalty_y,group_y=group_y,beta=beta,
                         max_iterations=max_iterations,tolerance=tolerance,
                         cross_validate=cross_validate,num_folds=num_folds,parallel_CV=parallel_CV,num_cores=num_cores,
                         comp
    )

    alpha <- out$alpha
    beta <- out$beta
    xi <- out$xi
    eta <- out$eta
    # make predictions (eta, xi) orthogonal to all other predictions
    if (comp > 1) {
      project_out <- function(mat, vec) {
        vec <- vec - mat %*% solve(t(mat) %*% mat, t(mat) %*% vec)
        vec/sqrt(sum(vec^2))
      }
      xi <- project_out(U_comps, xi)
      eta <- project_out(V_comps, eta)
    
      # get sparse wrt Xc,Yc if comp > 1
      #do enet or gglasso--use one se rule and extract best lambda and best coefs
      if (penalization_x=='none') {
        alpha <- out$alpha
      } else if (penalization_x=='enet') {
        refit_out_x <- glmnet::cv.glmnet(Xc, eta, intercept=FALSE)
        
        lam_x <- refit_out_x$lambda.1se
        alpha <- as.matrix(coef(refit_out_x, s = 'lambda.1se')[-1])
      } else if (penalization_x=='glasso') {
        refit_out_x <- gglasso::cv.gglasso(Xc, eta, group=group_x)
        
        lam_x <- refit_out_x$lambda.1se
        alpha <- as.matrix(coef(refit_out_x, s = 'lambda.1se')[-1])
      }
      
      xi <- Xc %*% alpha
      alpha <- alpha/sqrt(sum(xi^2)) ## changes to normalize alpha
      xi <- xi/sqrt(sum(xi^2)) # same as Xc %*% alpha
      
      if (sum(is.na(alpha))>0) {
        list_comps[[comp]] <- list('intial cv'=out, 
                                   'refit x on original coords'=refit_out_x, 
                                   'refit x lambda'=lam_x)
        warning(paste0('Component ', comp, ' was set to zero.'))
        break # out of comp loop
      }

      if (penalization_x == 'none') {
        beta <- out$beta
      } else if (penalization_x == 'enet') {
        refit_out_y <- glmnet::cv.glmnet(Yc, xi, intercept=FALSE)
        
        lam_y <- refit_out_y$lambda.1se
        beta <- as.matrix(coef(refit_out_y, s = 'lambda.1se')[-1])
      } else if (penalization_x == 'glasso') {
        refit_out_y <- gglasso::cv.gglasso(Yc, xi, group=group_y)
        
        lam_y <- refit_out_y$lambda.1se
        beta <- as.matrix(coef(refit_out_y, s = 'lambda.1se')[-1])
      }
      
      eta <- Yc %*% beta
      beta <- beta/sqrt(sum(eta^2)) ## changes to normalize alpha
      eta <- eta/sqrt(sum(eta^2))
      
      if (sum(is.na(beta))>0) {
        list_comps[[comp]] <- list('intial cv'=out, 
                                   'refit x on original coords'=refit_out_x, 'refit y on original coords'=refit_out_y, 
                                   'refit x lambda'=lam_x, 'refit y lambda'=lam_y)
        warning(paste0('Component ', comp, ' was set to zero.'))
        break # out of comp loop
      }
    } # end if comp > 1
    
    # update A/B and U/V_components
    A <- cbind(A, alpha)
    B <- cbind(B, beta)
    U_comps <- cbind(U_comps, xi)
    V_comps <- cbind(V_comps, eta)
    
    # deflate
    Xdef <- Xdef - xi %*% (t(xi) %*% Xdef)
    Ydef <- Ydef - eta %*% (t(eta) %*% Ydef)
    
    # update list with all the info
    if (comp > 1) {
      list_comps[[comp]] <- list('intial cv'=out, 
                               'refit x on original coords'=refit_out_x, 'refit y on original coords'=refit_out_y, 
                               'refit x lambda'=lam_x, 'refit y lambda'=lam_y)
    } else {
      list_comps[[comp]] <- list('intial cv'=out)
    }
  } # end comp loop
  
  # scale back to original scale
  nmlize <- function(a) a/sqrt(sum(a^2))

  sb_x <- function(alpha) {
    nmlize(alpha/attributes(scale(X))$'scaled:scale')
  }
  sb_y <- function(beta) {
    nmlize(beta/attributes(scale(Y))$'scaled:scale')
  }
  A <- apply(A, 2, sb_x)
  B <- apply(B, 2, sb_y)
  
  if (ncol(A)<num_components) { # just making sure A and B have the expected size
    A <- cbind(A, matrix(NA, nrow=nrow(A), ncol=num_components-ncol(A)))
    B <- cbind(B, matrix(NA, nrow=nrow(B), ncol=num_components-ncol(B)))
  }
  
  # return everything
  return(list('alpha'=A,'beta'=B,'all_other'=list_comps))
}

#### cca algorithm to do the optimization ####
cca_algorithm <- function(X, Y, 
                          alpha=NA, beta=NA, 
                          penalization_x='none', penalization_y='none', 
                          overallpenalty_x, overallpenalty_y, 
                          relpenalty_x, relpenalty_y, 
                          group_x, group_y, 
                          tolerance=10^-20, max_iterations=100
                          ) {
  if(any(is.na(alpha))) alpha <- matrix(rep(1, ncol(X))/sqrt(sum((X %*% rep(1, ncol(X)))^2)), ncol=1)
  if(any(is.na(beta))) beta <- matrix(rep(1, ncol(Y)) / sqrt(sum((Y %*% rep(1, ncol(Y)))^2)), ncol=1)
  
  coef_diffs <- c(Inf)
  iteration_count <- 1
  converged <- FALSE
  
  #xi <- scale(X %*% alpha)
  xi <- X %*% alpha
  alpha <- alpha/sqrt(sum(xi^2))
  xi <- xi/sqrt(sum(xi^2))
  
  #eta <- scale(Y %*% beta)
  eta <- Y %*% beta
  beta <- beta/sqrt(sum(eta^2))
  eta <- eta/sqrt(sum(eta^2))
  
  while(!converged && iteration_count < max_iterations) { # loop coordinate descent until stop conditions
    #print(paste0('iteration count: ', iteration_count))
    
    alpha_new <- switch(penalization_x,
                     "none" = {
                       solve(t(X) %*% X, t(X) %*% eta)
                     },
                     "enet" = {
                       as.matrix(get_enet(X, eta, overallpenalty_x, relpenalty_x))
                     },
                     "glasso" = {
                       as.matrix(get_glasso(X, eta, group_x, overallpenalty_x))
                     }
    )
    #xi <- scale(X %*% alpha_new)
    
    if (sum(abs(alpha_new))==0) { # this happens eg when alpha_new = rep(0, p)
      xi <- matrix(rep(0, length(xi)), ncol=1)
      eta <- matrix(rep(0, length(eta)), ncol=1)
      alpha <- matrix(rep(0, length(alpha)), ncol=1)
      beta <- matrix(rep(0, length(beta)), ncol=1)
      
      iteration_count <- iteration_count + 0.5
      break
    }
    
    beta_new <- switch(penalization_y,
                     "none" = {
                       solve(t(Y) %*% Y, t(Y) %*% xi)
                     },
                     "enet" = {
                       as.matrix(get_enet(Y, xi, overallpenalty_y, relpenalty_y))
                     },
                     "glasso" = {
                       as.matrix(get_glasso(Y, xi, group_y, overallpenalty_y))
                     }
    )
    #eta <- scale(Y %*% beta_new)
    
    if (sum(abs(beta_new))==0) {
      xi <- matrix(rep(0, length(xi)), ncol=1)
      eta <- matrix(rep(0, length(eta)), ncol=1)
      alpha <- matrix(rep(0, length(alpha)), ncol=1)
      beta <- matrix(rep(0, length(beta)), ncol=1)
      
      iteration_count <- iteration_count + 1
      break
    }
    
    # update stop conditions
    alpha <- alpha_new
    beta <- beta_new
    
    iteration_count <- iteration_count + 1
    coef_diffs[iteration_count] <- sqrt(sum((alpha - alpha_new)^2, (beta - beta_new)^2)/(length(alpha) + length(beta)))
    
    if (abs(coef_diffs[iteration_count] - coef_diffs[iteration_count-1]) < tolerance){
      converged <- TRUE
    }
    
    #xi <- scale(X %*% alpha)
    xi <- X %*% alpha
    alpha <- alpha/sqrt(sum(xi^2))
    xi <- xi/sqrt(sum(xi^2))
    
    #eta <- scale(Y %*% beta)
    eta <- Y %*% beta
    beta <- beta/sqrt(sum(eta^2))
    eta <- eta/sqrt(sum(eta^2))
    
  }# end of main loop for coord desc
  #print(paste0('number of iterations: ', iteration_count))
  
  rownames(alpha) <- colnames(X)
  rownames(beta) <- colnames(Y)
  
  return(list('alpha'=alpha,'beta'=beta,'xi'=xi,'eta'=eta,'num_iterations'=iteration_count))
}


#### copy and pasted all other functions from original push, barely edited ####

#### FUNCTIONS FOR ELASTIC NET PENALIZATION ####
get_enet = function(X, y, lambda, alpha){
  as.numeric(glmnet::glmnet(X, y, lambda = lambda, alpha = alpha, intercept=FALSE)$beta) # removed intercept
}
get_enet_grid <- function(mat, relpenalty, nlambda, lam_eps) {
  overallpenalty_m <- max(abs(mat))/relpenalty # actual max
  overallpenalty_m <- 1.5*overallpenalty_m # but scaling up could be helpful
  #print(overallpenalty_m)
  exp(seq(log(overallpenalty_m), log(lam_eps*overallpenalty_m), length.out = nlambda))
}

#### FUNCTION FOR GROUP LASSO PENALIZATION ####
get_glasso <- function(X, y, group, lambda) {
  beta <- c(gglasso::gglasso(X, y, group = group, lambda = lambda, intercept=FALSE)$beta) #removed intercept
  names(beta) <- colnames(X)
  return(beta)
}
get_glasso_grid <- function(mat, group, nlambda, lam_eps) {
  vec <- vector(length = max(group)) # assumes groups are numbered 1, 2, ..., g
  for (i in 1:length(vec)) {
    vec[i] <- sqrt(sum((mat[group==i])^2))
  }
  
  overallpenalty_m <- max(vec) # actual max
  overallpenalty_m <- 1.5*overallpenalty_m # but scaling up could be helpful
  exp(seq(log(overallpenalty_m), log(lam_eps*overallpenalty_m), length.out = nlambda))
}

#### FUNCTION TO CROSS VALIDATE PENALTY PARAMETERS ####
get_cross_validated_penalty_parameters <- function(X, Y,
                                                   penalization_x, penalization_y,
                                                   overallpenalty_x, overallpenalty_y,
                                                   relpenalty_x, relpenalty_y,
                                                   group_x, group_y,
                                                   num_folds, num_cores, parallel_CV,
                                                   max_iterations, tolerance
){
  shuffled <- get_split_sets(X = X, Y = Y, nr_subsets = num_folds)
  
  X.sampled <- shuffled$X.sampled
  Y.sampled <- shuffled$Y.sampled
  label <- shuffled$labels
  
  if (parallel_CV){
    print('have not yet updated parallel CV with new penalty params')
    cv_results <- get_parallel_cv(X = X.sampled, Y = Y.sampled,
                                  penalization_x = penalization_x, penalization_y = penalization_y,
                                  lasso_penalty_x = lasso_penalty_x, lasso_penalty_y = lasso_penalty_y,
                                  grp_penalty_x = grp_penalty_x, grp_penalty_y = grp_penalty_y,
                                  ridge_penalty_x = ridge_penalty_x, ridge_penalty_y = ridge_penalty_y,
                                  group_x = group_x, group_y = group_y,
                                  max_iterations = max_iterations, tolerance = tolerance,
                                  label = label, parallel_CV = parallel_CV, nr_cores = num_cores)
  } else {
    cv_results <- get_non_parallel_cv(X = X.sampled, Y = Y.sampled,
                                      penalization_x = penalization_x, penalization_y = penalization_y,
                                      overallpenalty_x = overallpenalty_x, overallpenalty_y = overallpenalty_y,
                                      relpenalty_x = relpenalty_x, relpenalty_y = relpenalty_y,
                                      group_x = group_x, group_y = group_y,
                                      max_iterations = max_iterations, tolerance = tolerance,
                                      label = label)
  }
  print("############Exited cross validation")
  
  # CV error minimizer
  #a = cv_results$mean_abs_cors[,'corr mean'] #BB: changed from 3 to 7 to 5
  #best_values <- cv_results$mean_abs_cors[which.max(a),]
  
  # CV error one-se minimizer
  mac <- cv_results$mean_abs_cors
  max_corr_ind <- which.max(mac[,'corr mean'])
  within_ose_mac <- mac[mac[,'corr mean'] > mac[max_corr_ind,'corr mean']-mac[max_corr_ind,'corr se'], ]
  best_ind <- which.max(within_ose_mac[,'sparsity']) # should never have length > 1
  best_values <- within_ose_mac[best_ind,]
  
  list(
    abs_cors = cv_results$abs_cors,
    mean_abs_cors = cv_results$mean_abs_cors,
    stime = cv_results$stime,
    num_iterations = cv_results$num_iterations,
    best_rel_y = best_values[1],
    best_overall_y = best_values[2],
    best_rel_x = best_values[3],
    best_overall_x = best_values[4]
  )
}


#### FUNCTION TO SPLIT UP DATA MATRIX INTO n SUBSETS ####
get_split_sets <- function(X, Y, nr_subsets){
  n <- dim(X)[1]
  
  #re-sample data rows
  splitting_dimensions <- sample(1:n,n)
  
  X.sampled <- X[splitting_dimensions,]
  Y.sampled <- Y[splitting_dimensions,]
  
  #calculate how many rows are left over
  leftover = n %% nr_subsets
  
  rep_eat = (n-leftover)/nr_subsets
  
  #repeat sequence nr_subsets
  labels = rep(1:nr_subsets, each=rep_eat)
  
  if(leftover!=0) {
    labels = c(labels, (1:nr_subsets)[1:leftover])
  }
  
  list(X.sampled = X.sampled, Y.sampled = Y.sampled, labels = labels)
}


#### CROSS VALIDATION FOR FINDING OPTIMAL PENALIZATION PARAMS ####
get_non_parallel_cv <- function(X, Y,
                                penalization_x, penalization_y,
                                overallpenalty_x, overallpenalty_y,
                                relpenalty_x, relpenalty_y,
                                group_x, group_y,
                                max_iterations, tolerance,
                                label){
  num_folds <- length(unique(label))
  abs_cors <- c()
  num_iterations <- c()
  sparsity <- c()
  kth_fold <- 0
  
  ii <- 0
  
  #overallpenalty_y <- c(.02, .06)
  #overallpenalty_x <- c(.02, .06)
  
  #l1x=1;l2x=1;l1y=1;l2y=1;
  stime <- system.time({
      for (l1x in 1:length(overallpenalty_x)){
        for (l2x in 1:length(relpenalty_x)){
          for (l1y in 1:length(overallpenalty_y)){
              for (l2y in 1:length(relpenalty_y)){
                ii <- ii+1 
                print(paste("#####Doing tuning param combo #", ii, '')) 
                
                sub_abs_cor <- c()
                sub_num_iterations <- c()
                sub_sparsity <- c()
                
                alpha.old <- NA
                beta.old <- NA
                for (i in 1:num_folds){
                  #print(paste0("doing fold number", i))
                  kth_fold <- kth_fold + 1
                  if (sum(abs(alpha.old),na.rm=TRUE)==0) alpha.old <- NA # if one fold killed coefs, next fold would give error
                  if (sum(abs(beta.old),na.rm=TRUE)==0) beta.old <- NA
                  
                  X.train <- X[label!=i,]
                  X.test <- X[label==i,]
                  Y.train <- Y[label!=i,]
                  Y.test <- Y[label==i,]
                  
                  sub_results <- cca_algorithm(X = X.train, Y = Y.train,
                                               alpha = alpha.old, beta = beta.old,
                                               penalization_x = penalization_x, penalization_y = penalization_y,
                                               overallpenalty_x = overallpenalty_x[l1x], overallpenalty_y = overallpenalty_y[l1y],
                                               relpenalty_x = relpenalty_x[l2x], relpenalty_y = relpenalty_y[l2y],
                                               group_x = group_x, group_y = group_y,
                                               tolerance = tolerance, max_iterations = max_iterations)
                  
                  alpha.old <- sub_results$alpha
                  beta.old <- sub_results$beta
                  
                  if(i==1) {
                    #print(c(overallpenalty_x[l1x], overallpenalty_y[l1y]))
                    #print(alpha.old)
                    #print(beta.old)
                  }

                  xi.test <- scale(X.test) %*% alpha.old
                  
                  #Divide with dim(Y.train)[2], exclude NA's from sum
                  sub_abs_cor[i] <- ifelse(sum(abs(xi.test))==0, 0, sum(abs(cor(xi.test,Y.test)),na.rm = TRUE)) #changed to be 0 if all NA
                  sub_num_iterations[i] <- sub_results$num_iterations
                  sub_sparsity[i] <- mean(c(mean(alpha.old==0), mean(beta.old==0)))
                } # end of loop over folds
                abs_cors <- cbind(abs_cors, sub_abs_cor) #matrix (nrows=num_folds, ncols=num_tuning_param_combos)
                num_iterations <- cbind(num_iterations, sub_num_iterations)
                sparsity <- cbind(sparsity, sub_sparsity)
              }
            }
          }
    } # end of tuning params loop
  })[3]
  
  labels_tuning <- as.matrix(expand.grid(relpenalty_y, overallpenalty_y, relpenalty_x, overallpenalty_x)) #######
  
  mean_abs_cors <- c()
  for (i in 1:nrow(labels_tuning)){
    v <- abs_cors[,i]
    v <- v[!is.na(v)]
    mean_abs_cors <- rbind(mean_abs_cors, #BB: this is the existing rows
                           c(labels_tuning[i,], mean(v), sd(v)/sqrt(length(v)), mean(sparsity[,i]))) #included SE too
  }
  rownames(mean_abs_cors) <- NULL
  colnames(mean_abs_cors) <- c("Y relative penalty", "Y overall penalty",
                               "X relative penalty", "X overall penalty",
                               "corr mean", "corr se", "sparsity")
  
  ## Return section ##
  list(abs_cors = abs_cors,
       mean_abs_cors = mean_abs_cors,
       stime = stime,
       num_iterations = num_iterations)
}

#### PARALLEL N-FOLD CROSS VALIDATION ####
get_parallel_cv <- function(X,
                            Y, ###
                            penalization_x, ###
                            lasso_penalty_x,
                            grp_penalty_x,
                            ridge_penalty_x,
                            group_x,
                            penalization_y, ###
                            lasso_penalty_y,
                            grp_penalty_y,
                            ridge_penalty_y,
                            group_y,
                            label, ###
                            max_iterations,
                            tolerance,
                            parallel_CV,
                            nr_cores){
  
  stop('Parallel has not been updated yet')
  nr_subsets <- length(unique(label))
  
  ## PARALLELIZATION ##
  cl <- parallel::makeForkCluster(nr_cores)
  doParallel::registerDoParallel(cl)
  
  #Measure time
  stime <- system.time({
    x <- foreach(grp_pen_x = grp_penalty_x, .export= c("sCCA", "get_glasso", "get_enet"), .combine = c) %dopar% { 
      abs_cors <- c()
      num_iterations  <- c()
      
      #for (l2x in 1:length(grp_penalty_x)){
      for (l1x in 1:length(lasso_penalty_x)){
        for (l3x in 1:length(ridge_penalty_x)){
          for (l2y in 1:length(grp_penalty_y)){
            for (l1y in 1:length(lasso_penalty_y)){
              for (l3y in 1:length(ridge_penalty_y)){
                
                sub_abs_cor <- c()
                sub_num_iterations <- c() 
                
                ALPHA.old <- c() 
                BETA.old <- c()
                
                for (i in 1:nr_subsets){
                  X.train <- X[label!=i,]
                  X.test <- X[label==i,]
                  
                  Y.train <- Y[label!=i,]
                  Y.test <- Y[label==i,]

                  sub_results <- sCCA(X = X.train,
                                      Y = Y.train,
                                      penalization_x = penalization_x, #######
                                      lasso_penalty_x = lasso_penalty_x[l1x],
                                      grp_penalty_x = grp_pen_x,
                                      ridge_penalty_x = ridge_penalty_x[l3x],
                                      group_x = group_x,
                                      ALPHA = ALPHA.old,
                                      penalization_y = penalization_y, #######
                                      lasso_penalty_y = lasso_penalty_y[l1y],
                                      grp_penalty_y = grp_penalty_y[l2y],
                                      ridge_penalty_y = ridge_penalty_y[l3y],
                                      group_y = group_y,
                                      BETA = BETA.old,
                                      max_iterations = max_iterations,
                                      tolerance = tolerance,
                                      cross_validate = FALSE)
                  
                  ALPHA.old <- sub_results$ALPHA 
                  BETA.old <- sub_results$BETA 
                  
                  XI.test = scale(X.test) %*% ALPHA.old
                  
                  sub_num_iterations[[i]] <- sub_results$nr_iterations ##NB this is an atomic vector (ult. length = nr_subsets)
                  sub_abs_cor[[i]] <- sum(abs(cor(XI.test,Y.test)), na.rm = TRUE) ##BB: why sum here?
                } #End of subset for loop
                
                abs_cors <- cbind(abs_cors, sub_abs_cor) ##BB: NB this is a matrix (ult. nr_subsets x (#tuningparamsbesidesgrpx))
                num_iterations <- cbind(num_iterations, sub_num_iterations)
              } #End of all other params loop
            }
          }
        }
      }
      list(abs_cors, num_iterations)
    } #End of lambda for loop
  })[3] #End of measure time
  
  #Stop cluster
  stopCluster(cl) 
  
  #mine out the iterations from vector x 
  abs_cors <- c()
  num_iterations <- c()
  j <- 0
  for (i in 1:length(grp_penalty_x)) {
    j <- j+1
    abs_cors <- cbind(abs_cors, x[[j]])
    j <- j+1
    num_iterations <- cbind(num_iterations, x[[j]])
  }
  
  #Arrange like in get_non_parallel 
  labels_tuning <- expand.grid(ridge_penalty_y, lasso_penalty_y, grp_penalty_y, 
                               ridge_penalty_x, lasso_penalty_x, grp_penalty_x) 
  
  mean_abs_cors <- c()
  for (i in 1:nrow(labels_tuning)){
    v <- abs_cors[,i]
    v <- v[!is.na(v)]
    mean_abs_cors <- rbind(mean_abs_cors,
                           c(labels_tuning[i,], mean(v), sd(v)/sqrt(length(v)))) 
  }
  
  rownames(mean_abs_cors)   <- NULL
  colnames(mean_abs_cors)   <- c("Y Ridge Penalty", "Y Lasso Penalty", "Y Group Penalty",
                                 "X Ridge Penalty", "X Lasso Penalty", "X Group Penalty",
                                 "mean Cor over CVs", "SE Cor over CVs")
  
  list(abs_cors = abs_cors, mean_abs_cors = mean_abs_cors, stime = stime)
}
