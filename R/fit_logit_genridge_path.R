

#' Ridge solution path for logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector
#' @param x A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real vector) Vector of regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param use_nn (boolean) Use nearest neighbors search as safety check on PQ bounds optimization in later iterations.
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
#' 
fit_logit_ridge_path <- function(y, X, type=c('NR','BL','PG','PQ'), 
                                 beta_start=NULL, lambda=NULL, eps=1e-10, 
                                 gamma=1.0, intercept=FALSE, phi=0.9, maxiter=1000L, 
                                 abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                                 verbose=FALSE, freq=10L, 
                                 solver=c("auto", "chol", "smw", "admm"),
                                 ctr_admm=set_ctr_admm()){
  
  # Check the bound type and QP method
  type <- match.arg(type)
  solver <- match.arg(solver)
  
  # Set the lambda vector
  if (is.null(lambda)) {
    lambda <- 10^seq(from=-4, to=+3, by=0.25)
  } else {
    lambda <- sort(lambda, decreasing=FALSE)
  }
  
  # Set the lambda dimension
  n <- nrow(X)
  p <- ncol(X)
  K <- length(lambda)
  
  # Set the initial beta
  if (is.null(beta_start)) {
    beta_start <- rep(0, p)
    if (intercept) {
      beta_start[1] <- qlogis(mean(y))
    }
  }
  
  # Set the beta path
  out <- list()
  out$lambdas <- lambda
  out$alpha   <- .0
  out$beta    <- matrix(NA, nrow=p, ncol=K)
  out$loglik  <- rep(NA, times=K)
  out$niter   <- rep(NA, times=K)
  out$exetime <- rep(NA, times=K)
  out$dev     <- rep(NA, times=K)
  out$edf     <- rep(NA, times=K)
  out$reml    <- rep(NA, times=K)
  out$gcv     <- rep(NA, times=K)
  out$aic     <- rep(NA, times=K)
  out$bic     <- rep(NA, times=K)
  
  # Progress bar
  if (!verbose) cat("\n|")
  
  # Solution path
  for (k in K:1) {
    if (!verbose) cat(".")
    
    # Set the starting time
    timek <- proc.time()
    
    # Fit the current model with warm-start initialization
    fitk <- fit_logit_ridge(y=y, X=X, type=type, beta_start=beta_start, lambda=lambda[k], 
                            eps=eps, phi=phi, intercept=intercept, maxiter=maxiter, 
                            abstol=abstol, reltol=reltol, etatol=etatol, 
                            verbose=verbose, freq=freq, solver=solver,
                            ctr_admm=ctr_admm)
    
    
    # Set the next initial estimate
    beta_start <- fitk$beta
    
    # Compute the sufficient statistics for the edf and the logdet
    eta <- as.vector(X %*% fitk$beta)
    pr <- 1/(1+exp(-eta))
    wts <- pmax(1e-8, pr*(1-pr))
    pL2 <- rep(lambda[k], p)
    if (intercept) pL2[1] <- eps
    
    # GoF and complexity measures
    if (p>n) {
      # cat("SMW")
      XPXt <- X %*% (t(X)/pL2)
      R <- chol(diag(1/wts) + XPXt)
      edf <- smw_chol_edf(R, XPXt, wts, pL2)
      logdet <- smw_chol_logdet(R, wts, pL2)
    } else {
      # cat("Chol")
      XtWX <- crossprod(X, wts*X)
      R <- chol(XtWX + diag(pL2))
      edf <- chol_edf(R, XtWX, exact=FALSE, rank=10, nsample=100)
      logdet <- chol_logdet(R) - p*log(lambda[k])
    }
    
    reml <- fitk$loglik - 0.5*logdet
    dev <- mean(binomial()$dev.resid(y, pr, 1), na.rm=TRUE)
    gcv <- dev / (1 - gamma*edf/n)
    aic <- dev + (edf/n) * 2
    bic <- dev + (edf/n) * log(n)
    
    # Current execution time
    timek <- (proc.time() - timek)[3]
    
    # Current result storage
    out$beta[,k]   <- fitk$beta
    out$loglik[k]  <- fitk$loglik
    out$niter[k]   <- fitk$niter
    out$dev[k]     <- dev
    out$edf[k]     <- edf
    out$logdet[k]  <- logdet
    out$gcv[k]     <- gcv
    out$reml[k]    <- reml
    out$aic[k]     <- aic
    out$bic[k]     <- bic
    out$exetime[k] <- timek
  }
  if (!verbose) cat("|\n")
  
  # Return the estimated path
  return(out)
}




#' Ridge solution path for logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector
#' @param x A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real vector) Vector of regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param use_nn (boolean) Use nearest neighbors search as safety check on PQ bounds optimization in later iterations.
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
#' 
fit_logit_spridge_path <- function(y, X, D, type=c('NR','BL','PG','PQ'), 
                                   beta_start=NULL, lambda=NULL, eps=1e-10, 
                                   gamma=1.0, intercept=FALSE, phi=0.9, maxiter=1000L, 
                                   abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                                   verbose=FALSE, freq=10L, use_nn=FALSE, 
                                   method=c("prjg", "dual", "admm"), 
                                   ctr_admm=set_ctr_admm(), 
                                   ctr_prjg=set_ctr_prjg(),
                                   ctr_dual=set_ctr_dual()){
  
  # Check the bound type and QP method
  type <- match.arg(type)
  method <- match.arg(method)
  
  # Set the lambda vector
  if (is.null(lambda)) {
    lambda <- 10^seq(from=-4, to=+3, by=0.25)
  } else {
    lambda <- sort(lambda, decreasing=FALSE)
  }
  
  # Set the lambda dimension
  n <- nrow(X)
  p <- ncol(X)
  K <- length(lambda)
  
  # Set the model matrices
  X <- Matrix::Matrix(X, sparse=TRUE)
  D <- Matrix::Matrix(D, sparse=TRUE)
  P <- Matrix::crossprod(D)
  
  # Set the initial beta
  if (is.null(beta_start)) {
    beta_start <- rep(0, p)
  }
  
  # Set the beta path
  out <- list()
  out$lambdas <- lambda
  out$alpha   <- .0
  out$beta    <- matrix(NA, nrow=p, ncol=K)
  out$niter   <- rep(NA, times=K)
  out$exetime <- rep(NA, times=K)
  out$loglik  <- rep(NA, times=K)
  out$dev     <- rep(NA, times=K)
  out$edf     <- rep(NA, times=K)
  out$logdet  <- rep(NA, times=K)
  out$gcv     <- rep(NA, times=K)
  out$reml    <- rep(NA, times=K)
  out$aic     <- rep(NA, times=K)
  out$bic     <- rep(NA, times=K)
  
  # Progress bar
  if (!verbose) cat("\n|")
  
  # Solution path
  for (k in K:1) {
    if (!verbose) cat(".")
    
    # Set the starting time
    timek <- proc.time()
    
    # Fit the current model with warm-start initialization
    fitk <- fit_logit_spridge(y=y, X=X, D=D, type=type, beta_start=beta_start, 
                              lambda=lambda[k], eps=eps, intercept=intercept, phi=phi,
                              maxiter=maxiter, abstol=abstol, reltol=reltol, etatol=etatol, 
                              verbose=FALSE, freq=freq, use_nn=use_nn, method=method,
                              ctr_admm=ctr_admm, ctr_prjg=ctr_prjg, ctr_dual=ctr_dual)
    
    # Set the next initial estimate
    beta_start <- fitk$beta
    
    # Compute the sufficient statistics for the edf and the logdet
    eta <- as.vector(X %*% fitk$beta)
    pr <- 1 / (1+exp(-eta))
    wts <- pmax(1e-8, pr*(1-pr))
    XtWX <- Matrix::crossprod(sqrt(wts)*X)
    R <- Matrix::Cholesky(XtWX+lambda[k]*P)
    
    # GoF and complexity measures
    dev <- mean(binomial()$dev.resid(y, pr, 1), na.rm=TRUE)
    edf <- sp_chol_edf(R, XtWX, exact=FALSE, rank=10, nsample=100)
    logdet <- sp_chol_logdet(R) - p*log(lambda[k])
    reml <- fitk$loglik - 0.5*logdet
    gcv <- dev / (1 - gamma * edf/n)
    aic <- dev + (edf/n) * 2 
    bic <- dev + (edf/n) * log(n)
    
    # Current execution time
    timek <- (proc.time() - timek)[3]
    
    # Current result storage
    out$beta[,k]   <- fitk$beta
    out$loglik[k]  <- fitk$loglik
    out$niter[k]   <- fitk$niter
    out$dev[k]     <- dev
    out$edf[k]     <- edf
    out$logdet[k]  <- logdet
    out$gcv[k]     <- gcv
    out$reml[k]    <- reml
    out$aic[k]     <- aic
    out$bic[k]     <- bic
    out$exetime[k] <- timek
  }
  if (!verbose) cat("|\n")
  
  # Return the estimated path
  return(out)
}



#' Generalized lasso cross-validation for logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector
#' @param x A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real vector) Vector of regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param use_nn (boolean) Use nearest neighbors search as safety check on PQ bounds optimization in later iterations.
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
#' 
fit_logit_spridge_cv <- function(y, X, D, type=c('NR','BL','PG','PQ'), 
                                 nfold=5, seed=1234, beta_start=NULL, 
                                 lambda=NULL, gamma=1.0, phi=0.9, maxiter=1000L, 
                                 abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                                 verbose=FALSE, freq=10L, use_nn=FALSE,
                                 method=c("prjg", "dual", "admm"), 
                                 ctr_admm=set_ctr_admm(), 
                                 ctr_prjg=set_ctr_prjg(),
                                 ctr_dual=set_ctr_dual()) {
  
  # Set the seed for RNG
  set.seed(seed)
  
  # Split the data into nfold sub-groups
  n <- nrow(X)
  ridx <- permute::shuffle(n)
  fidx <- rep(1:nfold, length=n)
  folds <- split(ridx, fidx)
  
  # Cross-validation cycle
  df <- data.frame()
  for (k in 1:nfold) {
    # Report the current fold
    cat("fold:", k, "\n")
    
    # Split the train and test data
    train <- do.call(c, folds[-k])
    test <- folds[[k]]
    
    # Get the train and test dimension
    n_tr <- length(train)
    n_ts <- length(test)
    
    # Re-scale the lambdas wtr the training sample size
    lambdak <- lambda*(n_tr/n)
    
    # Fit the model on the training set
    fit <- fit_logit_spridge_path(y[train], X[train,], D, type=type, 
                                  beta_start=beta_start, lambda=lambdak, 
                                  gamma=gamma, phi=phi, maxiter=maxiter, 
                                  abstol=abstol, reltol=reltol, etatol=etatol, 
                                  verbose=verbose, freq=freq, use_nn=use_nn, method=method,
                                  ctr_admm=ctr_admm, ctr_prjg=ctr_prjg, ctr_dual=ctr_dual)
    
    # Compute the log-likelihood on the test set
    eta_ts <- as.matrix(X[test,] %*% fit$beta)
    pr_ts <- 1 / (1+exp(-eta_ts))
    loglik <- (n/n_ts) * colSums(y[test]*eta_ts - Rmpfr::log1pexp(eta_ts))
    
    # Get the iterations and execution times to convergence
    niter <- fit$niter
    exetime <- fit$exetime
    
    # Store all the  summary information
    dfk <- cbind(fold=k, lambda=lambda, loglik=loglik, niter=niter, exetime=exetime)
    df <- rbind(df, data.frame(dfk))
  }
  
  # Select the optimal lambda
  loglikcv <- aggregate(loglik ~ lambda, data=df, FUN=mean)$loglik
  lambda_best <- lambda[which.max(loglikcv)]
  
  # Refit the model with the selected lambda
  cat("fold: refit \n")
  ctr_admm$rho <- ctr_admm$rho * lambda_best
  fit <- fit_logit_spridge(y, X, D, type=type, beta_start=beta_start, 
                           lambda=lambda_best, phi=phi, maxiter=maxiter, 
                           abstol=abstol, reltol=reltol, etatol=etatol, 
                           verbose=FALSE, freq=freq, use_nn=use_nn, method=method,
                           ctr_admm=ctr_admm, ctr_prjg=ctr_prjg, ctr_dual=ctr_dual)
  
  # Append the cross-validation results to the fitted model
  fit$cv <- df
  
  # Return the cross-validation results 
  return(fit)
}






