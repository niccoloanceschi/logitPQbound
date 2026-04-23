

#' @title Ridge logistic regression via MM optimization
#'
#' @param y  A \verb{n} dimensional binary vector
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real) Regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param phi (real) inertia parameter for the active set detection
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param solver (string) Method to be used for the solution of the inner optimization.
#' @param ctr_prjg (list) Control parameters for the projected gradient algorithm
#' @param ctr_dual (list) Control parameters for the dual QP algorithm
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_ridge <- function(y, X, type=c('NR','BL','PG','PQ'), 
                            beta_start=NULL, lambda=NULL, eps=1e-10, 
                            intercept=FALSE, phi=0.9, maxiter=1000L, 
                            abstol=1e-4, reltol=1e-4, etatol=1e-1, 
                            verbose=FALSE, freq=10L,
                            solver=c("auto", "chol", "smw", "admm"),
                            ctr_admm=set_ctr_admm()){
  
  # Check the bound type and QP method
  type = match.arg(type)
  solver = match.arg(solver)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  # Set the linear system solver
  if (solver=="prjg")
    stop("PRJG-SOLVER not implemented yet.", call.=FALSE)
  if (solver=="dual")
    stop("DUAL-SOLVER not implemented yet.", call.=FALSE)
  if (solver=="dual")
    stop("SPARSE-SOLVER not implemented yet.", call.=FALSE)
  
  if (type=="BL" & solver=="admm")
    warning("ADMM-SOLVER is not available for the BL bound.", immediate.=TRUE, call.=FALSE)
  if (type=="PG" & solver=="admm")
    warning("ADMM-SOLVER is not available for the PG bound.", immediate.=TRUE, call.=FALSE)
  if (type=="NR" & solver=="admm")
    warning("ADMM-SOLVER is not available for the NR bound.", immediate.=TRUE, call.=FALSE)
  if (type!="PQ" & solver=="admm")
    solver <- "auto"
  
  if (solver=="auto")
    solver <- ifelse(p<=n, "chol", "smw")
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("RIDGE-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda, p)
  if (intercept) pL2[1] <- eps
  
  # Set the sufficient statistics
  XtWy <- crossprod(X, y-0.5)
  XPXt <- NULL
  cholQ <- NULL
  if (p>n) {
    XP <- t(t(X)/pL2)
    XPXt <- tcrossprod(X, XP)
  }
  
  # Set the log-posterior function
  logpost <- function(y, eta, beta, pL1, pL2){
    loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
    pridge <- 0.5*sum(pL2*beta^2)
    plasso <- sum(pL1*abs(beta))
    return(loglik - pridge - plasso)
  }
  
  # Bound initialization
  if(is.null(beta_start)){beta_start <- rep(0,p)}
  thr <- 1e-2
  beta <- beta_start
  eta <- as.vector(X %*% beta)
  sign_eta <- soft_sign(eta, thr)
  
  # Set the initial optimization state
  iter <- 1
  objval <- logpost(y, eta, beta, .0, pL2)
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_beta <- NaN
  diff_rel_eta <- NaN
  
  # Bound initialization
  eta_stable <- FALSE
  if (type=='PQ') {
    coeff_pq <- coeff_PQ(eta)
    w  <- c(coeff_pq[,1])
    nu <- c(coeff_pq[,2])
    s  <- sign_eta
  } else if (type=='PG') {
    w <- c(coeff_PG(eta))
  } else if (type=='NR') { 
    step <- 1.0
    w <- exp(-eta)/(1+exp(-eta))^2
  } else if (type=='BL') {
    w <- 0.25
    if (p> n) cholQ <- chol(XPXt + diag(1/w,n,n))
    if (p<=n) cholQ <- chol(w*crossprod(X) + diag(pL2,p,p))
  }

  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=7))
  colnames(trace) <- c("iter", "objval", "diff_abs_obj", "diff_rel_obj", 
                       "diff_beta", "diff_eta", "stable")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, 
                    diff_rel_beta, diff_rel_eta, eta_stable)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "beta"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "stable"), "\n",
        paste0(rep("-", 73), collapse=""), "\n",
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
    
    # Set the old parameters
    objval_old <- objval
    beta_old <- beta 
    eta_old <- eta
    sign_eta_old <- sign_eta
    
    # Relax the ADMM tolerance for the first iterations
    if (TRUE) {
      ctr_admm$objtol <- objtol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$abstol <- abstol_admm * (1+sqrt(p/n)*100*iter^-.75)
      ctr_admm$reltol <- reltol_admm * (1+sqrt(p/n)*100*iter^-.75)
    }
    
    # Switch to the stationary update when sign(eta) is stable
    if (eta_stable & solver=="admm") {
      solver <- ifelse(p<=n, "chol", "smw")
    }
    
    # Update beta
    beta <- switch(type,
      'NR' = update_ridge_NR(y, X, eta, w, pL2, XPXt, XtWy, solver),
      'BL' = update_ridge_BL(y, X, eta, pL2, cholQ, solver),
      'PG' = update_ridge_PG(y, X, eta, w, pL2, XtWy, XPXt, solver),
      'PQ' = update_ridge_PQ(y, X, eta, w, nu, s, pL2, beta, phi, XPXt, solver, ctr_admm))
    
    # Update bound
    eta <- as.vector(X %*% beta)
    sign_eta <- soft_sign(eta, thr)
    
    # Check convergence
    objval <- logpost(y, eta, beta, .0, pL2)
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_beta <- maxreldiff(beta, beta_old)
    diff_rel_eta <- maxreldiff(eta, eta_old)
    min_prod_eta <- min(sign_eta*sign_eta_old) # min(eta*eta_old)
    
    if (!eta_stable) {
      eta_stable <- (diff_rel_eta<etatol) && (min_prod_eta>=0)
    }
    
    if (type=='PQ'){         
      coeff_pq <- coeff_PQ(eta)
      w  <- c(coeff_pq[,1])
      nu <- c(coeff_pq[,2])
      s  <- (1-phi)*s + phi*sign_eta
      s  <- ifelse(sign_eta>thr, sign_eta, s)
    } else if (type=='PG'){ 
      w <- c(coeff_PG(eta))
    } else if (type=='BL'){
      w <- rep(0.25, n)
    } else if (type=='NR'){
      w <- exp(-eta)/(1+exp(-eta))^2
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
        gettextf("%10s", eta_stable), "\n",
        paste0(rep("-", 73), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("RIDGE-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  return(list(beta=beta, niter=iter, loglik=objval, weights=w, trace=trace[1:iter,]))
}


#' @title Ridge logistic regression via MM optimization
#'
#' @param y A \verb{n} dimensional binary vector
#' @param x A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param D A \verb{q x p} dimensional penalty matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real) Regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param method (string) Method to be used for the solution of the quadratic programming inner optimization. Must be one of 'dual', or 'admm'
#' @param ctr_admm (list) Control parameters for the ADMM algorithm
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_spridge <- function(y, X, D, type=c('NR','BL','PG','PQ'), 
                              beta_start=NULL, lambda=NULL, eps=1e-10, 
                              intercept=FALSE, phi=0.9, maxiter=1000L, 
                              abstol=1e-4, reltol=1e-4, etatol=1e-4, 
                              verbose=FALSE, freq=10L, 
                              ctr_admm=set_ctr_admm()){
  
  # Check the bound type and QP method
  type = match.arg(type)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("SP-RIDGE-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda, p)
  if (intercept) pL2[1] <- eps
  
  # Set the sufficient statistics
  X <- Matrix::Matrix(X, sparse=TRUE)
  D <- Matrix::Matrix(D, sparse=TRUE)
  P <- Matrix::crossprod(D)
  
  # Set the log-posterior function
  logpost <- function(beta, y, X, D, lambda, alpha){
    Xb <- as.vector(X %*% beta)
    Db <- as.vector(D %*% beta)
    loglik <- sum(y*Xb - Rmpfr::log1pexp(Xb))
    pridge <- lambda*(1-alpha)*0.5*sum(Db*Db)
    plasso <- lambda*alpha*sum(abs(Db))
    return(loglik - pridge - plasso)
  }
  
  # Bound initialization
  if(is.null(beta_start)){beta_start <- rep(0,p)}
  thr <- 1e-2
  beta <- as.vector(beta_start)
  eta <- as.vector(X %*% beta)
  sign_eta <- soft_sign(eta, thr)
  
  # Set the initial optimization state
  iter <- 1
  objval <- logpost(beta, y, X, D, lambda, .0)
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_beta <- NaN
  diff_rel_eta <- NaN
  
  # Set some control parameter
  thr <- 1e-2 # active set threshold
  step <- 0.99 # Newton step-size
  
  # Bound initialization
  eta_stable <- FALSE
  if (type=='PQ') {
    coeff_pq <- coeff_PQ(eta)
    w  <- c(coeff_pq[,1])
    nu <- c(coeff_pq[,2])
    s  <- sign_eta
  } else if (type=='PG') {
    Xtz <- Matrix::crossprod(X, y-0.5)
    w <- c(coeff_PG(eta))
  } else if (type=='NR') {      
    prob <- 1/(1+exp(-eta))
    w <- prob*(1-prob)
  } else if (type=='BL') {
    prob <- 1/(1+exp(-eta))
    w <- rep(0.25, times=n)
    Q <- Matrix::crossprod(rbind(sqrt(w)*X, sqrt(lambda)*D))
    R <- Matrix::Cholesky(Q, pivot=TRUE)
  }
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=7))
  colnames(trace) <- c("iter", "objval", "diff_abs_obj", "diff_rel_obj", 
                       "diff_beta", "diff_eta", "stable")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, 
                    diff_rel_beta, diff_rel_eta, eta_stable)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "beta"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "stable"), "\n",
        paste0(rep("-", 73), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/n),
        "\n", sep="")
  }
  
  # Iterative procedure
  for (iter in 2:maxiter) {
    
    # Set the old parameters
    objval_old <- objval
    beta_old <- beta 
    eta_old <- eta
    
    # Update beta
    if (type=="NR") {
      Q <- Matrix::crossprod(rbind(sqrt(w)*X, sqrt(lambda)*D))
      R <- Matrix::Cholesky(Q, pivot=TRUE)
      Xtz <- Matrix::crossprod(X, y-prob) - lambda * (P %*% beta)
      beta <- beta + step*sp_chol_solve(R, Xtz)
    }
    if (type=="BL") {
      Xtz <- Matrix::crossprod(X, y-prob) - lambda * (P %*% beta)
      beta <- beta + sp_chol_solve(R, Xtz)
    }
    if (type=="PG") {
      Q <- Matrix::crossprod(rbind(sqrt(w)*X, sqrt(lambda)*D))
      R <- Matrix::Cholesky(Q, pivot=TRUE)
      beta <- sp_chol_solve(R, Xtz)
    }
    if (type=="PQ") {
      Q <- Matrix::crossprod(rbind(sqrt(w)*X, sqrt(lambda)*D))
      R <- Matrix::Cholesky(Q, pivot=TRUE)
      Xtz <- Matrix::crossprod(X, y - 0.5 - nu*s)
      beta <- sp_chol_solve(R, Xtz)
    }
    
    # Update bound
    eta <- as.vector(X %*% beta)
    sign_eta <- soft_sign(eta, thr)
    
    # Check convergence
    objval <- logpost(beta, y, X, D, lambda, .0)
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_beta <- maxreldiff(beta, beta_old)
    diff_rel_eta <- maxreldiff(eta, eta_old)
    min_prod_eta <- min(eta*eta_old)
    eta_stable <- (diff_rel_eta<etatol) && (min_prod_eta>0)
    
    if (type=='PQ') {
      coeff_pq <- coeff_PQ(eta)
      w  <- c(coeff_pq[,1])
      nu <- c(coeff_pq[,2])
      s  <- (1-phi)*s + phi*sign_eta
      s  <- ifelse(sign_eta>thr, sign_eta, s)
    } else if (type=='PG'){ 
      w <- c(coeff_PG(eta))
    } else if (type=='BL'){
      prob <- 1/(1+exp(-eta))
      w <- rep(0.25, n)
    } else if (type=='NR'){ 
      prob <- 1/(1+exp(-eta))
      w <- prob*(1-prob)
    }
    
    # Store the current state
    trace[iter,] <- c(iter, objval, diff_abs_obj, diff_rel_obj, 
                      diff_rel_beta, diff_rel_eta, eta_stable)
    
    # Print the current optimization state
    if (verbose && iter%%freq==0) {
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/n),
          gettextf(" %10.3e", diff_abs_obj/n),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_beta),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf("%10s", eta_stable),
          "\n", sep="")
    }
    
    if ((diff_abs_obj<n*abstol) && (diff_rel_obj<reltol)) {
      break
    }
  }
  
  # Print the final optimization state
  if (verbose) {
    cat(gettextf("%6d", iter),
        gettextf("%12.5f", objval/n),
        gettextf("%11.3e", diff_abs_obj/n),
        gettextf("%11.3e", diff_rel_obj),
        gettextf("%11.3e", diff_rel_beta),
        gettextf("%11.3e", diff_rel_eta),
        gettextf("%10s", eta_stable), "\n",
        paste0(rep("-", 73), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("SP-RIDGE-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  return(list(beta=beta, niter=iter, loglik=objval, weights=w, trace=trace[1:iter,]))
}
