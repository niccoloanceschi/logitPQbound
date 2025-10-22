

#' @title Generalized lasso logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector.
#' @param X A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param D A \verb{m x p} dimensional penalty matrix.
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param lambda (real) Regularization parameter for elastic-net regression.
#' @param alpha (real) Mixing parameter regulating the balancing between lasso and ridge penalties.
#' @param eps (real) Small value added to the intercept term to avoid numerical issues.
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error.
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood.
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood.
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor.
#' @param verbose (boolean) Print the intermediate state of the optimization.
#' @param freq (int) How often print the optimization state.
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'.
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_genlasso <- function(y, X, D, type=c('NR','BL','PG','PQ'), beta_start=NULL, 
                               lambda=1.0, alpha=.99, eps=1e-10, intercept=FALSE,
                               maxiter=1000, abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                               verbose=FALSE, freq=10, ctr_admm=set_ctr_admm()){
  
  # Check the bound type and QP method
  type = match.arg(type)
  method = match.arg(method)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("GENLASSO-LOGIT: lambda must be a positive real value.", call.=FALSE)
  if ((alpha>1) || (alpha<0))
    stop("GENLASSO-LOGIT: alpha must be a value in the interval [0,1].", call.=FALSE)
  
  pL2 <- rep(lambda*(1-alpha), p)
  pL1 <- rep(lambda*alpha, m)
  
  if (alpha==0) pL1 <- eps
  if (intercept) pL2[1] <- eps*lambda*alpha
  
  # Parameter initialization
  if(is.null(beta_start)){
    beta_start <- init_genlasso_beta(y, X, D, pL1, pL2)
  }
  
  beta <- as.vector(beta_start)
  eta <- as.vector(X %*% beta)
  
  # Bound initialization
  eta_stable <- FALSE
  if (type=='PQ') {
    coeff_pq <- coeff_PQ(eta)
    w  <- as.vector(coeff_pq[,1])
    nu <- as.vector(coeff_pq[,2])
    u <- numeric(n+m)
    z <- c(nu * as.vector(X %*% beta), as.vector(D %*% beta))
  } else if (type=='PG') {
    w <- as.vector(coeff_PG(eta))
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  } else if (type=='NR') { 
    prob <- 1/(1+exp(-eta))
    w <- prob*(1-prob)
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  } else if (type=='BL') {
    prob <- 1/(1+exp(-eta))
    w <- 0.25
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  }
  
  # Set the log-posterior function
  logpost <- function(y, eta, beta, D, pL1, pL2){
    z <- abs(as.vector(D %*% beta))
    loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
    pridge <- sum(pL2*z^2)/2
    plasso <- sum(pL1*z)
    return(loglik - pridge - plasso)
  }
  
  # Set the initial optimization state
  iter <- 1
  objval <- logpost(y, eta, beta, D, pL1, pL2)
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_beta <- NaN
  diff_rel_eta <- NaN
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=7))
  colnames(trace) <- c("iter", "objval", "diff_abs_obj", "diff_rel_obj", "diff_beta", "diff_eta", "stable_eta")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, diff_rel_beta, diff_rel_eta, eta_stable)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "beta"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "stable"),
        gettextf("%10s", "niter"), "\n",
        paste0(rep("-", 83), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/n),
        "\n", sep="")
  }
  
  # Set the final tolerance values for the ADMM algorithm
  if (TRUE) {
    objtol_admm <- ctr_admm$objtol
    abstol_admm <- ctr_admm$abstol
    reltol_admm <- ctr_admm$reltol
  }
  
  # Iterative procedure
  for(iter in 2:maxiter){
    
    ## Relax the ADMM tolerance for the first iterations
    if (TRUE) {
      ctr_admm$objtol <- objtol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$abstol <- abstol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$reltol <- reltol_admm * (1+sqrt(p/n)*100*iter^-.75)
    }
    
    ## Update beta
    admmbeta <- switch(type,
      'NR' = update_genlasso_NR(beta, z, u, eta, prob, y, X, w, D, lambda, alpha, eps, ctr_admm),
      'BL' = update_genlasso_BL(beta, z, u, eta, prob, y, X, D, lambda, alpha, eps, ctr_admm),
      'PG' = update_genlasso_PG(beta, z, u, y, X, w, D, lambda, alpha, eps, ctr_admm),
      'PQ' = update_genlasso_PQ(beta, z, u, eta, eta_stable, y, X, w, nu, D, lambda, alpha, eps, ctr_admm)
    )
    
    beta_old <- beta
    niter <- admmbeta$niter
    beta <- admmbeta$beta
    u <- admmbeta$u
    z <- admmbeta$z
    
    # Update bound
    eta_old <- eta
    eta <- as.vector(X %*% beta)
    
    # Check convergence
    objval_old <- objval
    objval <- logpost(y, eta, beta, D, pL1, pL2)
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- abs(objval-objval_old) / (abs(objval_old)+1e-8)
    diff_rel_beta <- max(abs(beta-beta_old) / (abs(beta_old)+1e-8))
    diff_rel_eta <- max(abs(eta-eta_old) / (abs(eta_old)+1e-8))
    min_prod_eta <- min(eta*eta_old)
    eta_stable <- (diff_rel_eta<etatol) && (min_prod_eta>0)
    
    # Update weights
    if (type=='PQ'){         
      coeff_pq <- coeff_PQ(eta)
      w  <- as.vector(coeff_pq[,1])
      nu <- as.vector(coeff_pq[,2])
    } else if (type=='PG'){ 
      w <- as.vector(coeff_PG(eta))
    } else if (type=='BL'){
      prob <- 1/(1+exp(-eta))
      w <- 0.25
    } else if (type=='NR'){ 
      prob <- 1/(1+exp(-eta))
      w <- prob*(1-prob)
    }
    
    # Store the current state
    trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, 
                      diff_rel_beta, diff_rel_eta, eta_stable)
    
    # Print the current optimization state
    if(verbose && iter%%freq==0){
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/n),
          gettextf(" %10.3e", diff_abs_obj/n),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_beta),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf("%10s", eta_stable),
          gettextf("%10d", niter),
          "\n", sep="")
    }
    
    if((diff_abs_obj<n*abstol) && (diff_rel_obj<reltol)){
      break
    }
  }
  
  # Print the final optimization state
  if(verbose){
    cat(gettextf("%6d", iter),
        gettextf("%12.5f", objval/n),
        gettextf("%11.3e", diff_abs_obj/n),
        gettextf("%11.3e", diff_rel_obj),
        gettextf("%11.3e", diff_rel_beta),
        gettextf("%11.3e", diff_rel_eta),
        gettextf("%10s", eta_stable),
        gettextf("%10d", niter), "\n",
        paste0(rep("-", 83), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("GENLASSO-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  if (length(z)>m) {
    z <- z[-(1:n)]
    u <- u[-(1:n)]
  }
  
  return(list(type=type, alpha=alpha, lambda=lambda, beta=beta, z=z, u=u, w=w,
              niter=iter, loglik=objval, trace=trace[1:iter,]))
}




#' @title Generalized lasso logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector.
#' @param X A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param D A \verb{m x p} dimensional penalty matrix.
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param lambda (real) Regularization parameter for elastic-net regression.
#' @param alpha (real) Mixing parameter regulating the balancing between lasso and ridge penalties.
#' @param eps (real) Small value added to the intercept term to avoid numerical issues.
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error.
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood.
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood.
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor.
#' @param verbose (boolean) Print the intermediate state of the optimization.
#' @param freq (int) How often print the optimization state.
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'.
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_splasso <- function(y, X, D, type=c('NR','BL','PG','PQ'), beta_start=NULL, 
                              lambda=1.0, alpha=.99, eps=1e-10, intercept=FALSE, 
                              phi=0.9, approx=TRUE, maxiter=1000, 
                              abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                              verbose=FALSE, freq=10, ctr_admm=set_ctr_admm()){
  
  # Check the bound type and QP method
  type = match.arg(type)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("SP-LASSO-LOGIT: lambda must be a positive real value.", call.=FALSE)
  if ((alpha>1) || (alpha<0))
    stop("SP-LASSO-LOGIT: alpha must be a value in the interval [0,1].", call.=FALSE)
  
  pL2 <- rep(lambda*(1-alpha), p)
  pL1 <- rep(lambda*alpha, m)
  
  if (alpha==0) pL1 <- eps
  if (intercept) pL2[1] <- eps*lambda*alpha
  
  # Set the sufficient statistics
  X <- Matrix::Matrix(X, sparse=TRUE)
  D <- Matrix::Matrix(D, sparse=TRUE)
  
  # Parameter initialization
  if(is.null(beta_start)){
    beta_start <- init_splasso_beta(y, X, D, pL1, pL2)
  }
  
  beta <- as.vector(beta_start)
  eta <- as.vector(X %*% beta)
  
  # Bound initialization
  eta_stable <- FALSE
  if (type=='PQ') {
    coeff_pq <- coeff_PQ(eta)
    s <- sign(eta) * (abs(eta) > 1e-2)
    w  <- as.vector(coeff_pq[,1])
    nu <- as.vector(coeff_pq[,2])
    if (approx) {
      u <- numeric(m)
      z <- as.vector(D %*% beta)
    } else {
      u <- numeric(n+m)
      z <- c(nu*as.vector(X %*% beta), as.vector(D %*% beta))
    }
  } else if (type=='PG') {
    w <- as.vector(coeff_PG(eta))
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  } else if (type=='NR') { 
    prob <- 1/(1+exp(-eta))
    w <- prob*(1-prob)
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  } else if (type=='BL') {
    prob <- 1/(1+exp(-eta))
    w <- 0.25
    u <- numeric(m)
    z <- as.vector(D %*% beta)
  }
  
  # Set the log-posterior function
  logpost <- function(y, eta, beta, D, pL1, pL2){
    z <- abs(as.vector(D %*% beta))
    loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
    pridge <- sum(pL2*z^2)/2
    plasso <- sum(pL1*z)
    return(loglik - pridge - plasso)
  }
  
  # Set the initial optimization state
  iter <- 1
  objval <- logpost(y, eta, beta, D, pL1, pL2)
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_beta <- NaN
  diff_rel_eta <- NaN
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=7))
  colnames(trace) <- c("iter", "objval", "diff_abs_obj", "diff_rel_obj", "diff_beta", "diff_eta", "stable_eta")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, diff_rel_beta, diff_rel_eta, eta_stable)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "beta"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "stable"),
        gettextf("%10s", "niter"), "\n",
        paste0(rep("-", 83), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/n),
        "\n", sep="")
  }
  
  # Set the final tolerance values for the ADMM algorithm
  if (TRUE) {
    objtol_admm <- ctr_admm$objtol
    abstol_admm <- ctr_admm$abstol
    reltol_admm <- ctr_admm$reltol
  }
  
  # Iterative procedure
  for(iter in 2:maxiter){
    
    ## 
    if (TRUE) {
      ctr_admm$objtol <- objtol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$abstol <- abstol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$reltol <- reltol_admm * (1+sqrt(p/n)*100*iter^-.75)
    }
    
    ## Update beta
    admmbeta <- switch(type,
      'NR' = update_splasso_NR(beta, z, u, eta, prob, y, X, w, D, lambda, alpha, eps, ctr_admm),
      'BL' = update_splasso_BL(beta, z, u, eta, prob, y, X, D, lambda, alpha, eps, ctr_admm),
      'PG' = update_splasso_PG(beta, z, u, y, X, w, D, lambda, alpha, eps, ctr_admm),
      'PQ' = update_splasso_PQ(beta, z, u, s, eta, y, X, w, nu, D, lambda, alpha, eps, phi, approx, ctr_admm)
    )
    
    beta_old <- beta
    niter <- admmbeta$niter
    beta <- admmbeta$beta
    u <- admmbeta$u
    z <- admmbeta$z
    s <- admmbeta$s
    
    # Update bound
    eta_old <- eta
    eta <- as.vector(X %*% beta)
    
    # Check convergence
    objval_old <- objval
    objval <- logpost(y, eta, beta, D, pL1, pL2)
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_beta <- maxreldiff(beta, beta_old)
    diff_rel_eta <- maxreldiff(eta, eta_old)
    min_prod_eta <- min(eta*eta_old)
    eta_stable <- (diff_rel_eta<etatol) && (min_prod_eta>0)
    
    # Update weights
    if (type=='PQ'){         
      coeff_pq <- coeff_PQ(eta)
      w  <- as.vector(coeff_pq[,1])
      nu <- as.vector(coeff_pq[,2])
    } else if (type=='PG'){ 
      w <- as.vector(coeff_PG(eta))
    } else if (type=='BL'){
      prob <- 1/(1+exp(-eta))
      w <- 0.25
    } else if (type=='NR'){ 
      prob <- 1/(1+exp(-eta))
      w <- prob*(1-prob)
    }
    
    # Store the current state
    trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, 
                      diff_rel_beta, diff_rel_eta, eta_stable)
    
    # Print the current optimization state
    if(verbose && iter%%freq==0){
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/n),
          gettextf(" %10.3e", diff_abs_obj/n),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_beta),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf("%10s", eta_stable),
          gettextf("%10d", niter),
          "\n", sep="")
    }
    
    if((diff_abs_obj<n*abstol) && (diff_rel_obj<reltol)){
      break
    }
  }
  
  # Print the final optimization state
  if(verbose){
    cat(gettextf("%6d", iter),
        gettextf("%12.5f", objval/n),
        gettextf("%11.3e", diff_abs_obj/n),
        gettextf("%11.3e", diff_rel_obj),
        gettextf("%11.3e", diff_rel_beta),
        gettextf("%11.3e", diff_rel_eta),
        gettextf("%10s", eta_stable),
        gettextf("%10d", niter), "\n",
        paste0(rep("-", 83), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("SP-LASSO-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  if (length(z)>m) {
    z <- z[-(1:n)]
    u <- u[-(1:n)]
  }
  
  return(list(type=type, alpha=alpha, lambda=lambda, beta=beta, z=z, u=u, s=s,
              niter=iter, loglik=objval, trace=trace[1:iter,]))
}
