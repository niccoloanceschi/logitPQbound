

#' @title Ridge logistic regression via MM optimization
#'
#' @param y  A \verb{n} dimensional binary vector
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real) Regularization parameters for ridge regression
#' @param alpha (real) Mixing parameter regulating the balancing between lasso and ridge penalties.
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param phi (real) inertia parameter for the active set detection
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param etatol (real) Convergence threshold for the relative change in the linear predictor
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' @param use_nn (boolean) Use nearest neighbors search as safety check on PQ bounds optimization in later iterations.
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_enet_coord <- function(y, X, type=c('NR','BL','PG', 'PQ'), 
                                 beta_start=NULL, lambda=NULL, alpha=.99, eps=1e-10, 
                                 intercept=FALSE, phi=0.9, maxiter=1000L, 
                                 abstol=1e-4, reltol=1e-4, etatol=1e-1, 
                                 verbose=FALSE, freq=10L){
  
  # Check the bound type and QP method
  type = match.arg(type)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("ENET-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda*(1-alpha), p)
  pL1 <- rep(lambda*alpha, p)
  
  if (intercept) pL2[1] <- eps*(1-alpha)
  if (intercept) pL1[1] <- eps*alpha
  
  if (intercept) {
    X[, 1] <- 1
    X[,-1] <- scale(X[,-1], center=TRUE, scale=TRUE)
  }
  
  # Set the log-posterior function
  logpost <- function(y, eta, beta, pL1, pL2){
    loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
    pridge <- 0.5*sum(pL2*beta^2)
    plasso <- sum(pL1*abs(beta))
    return(loglik - pridge - plasso)
  }
  
  # Bound initialization
  if (is.null(beta_start)) {beta_start <- rep(0,p)}
  if (intercept) {beta_start[1] <- qlogis(mean(y))}
    
  thr <- 1e-2
  phi0 <- phi
  beta <- beta_start
  eta <- as.vector(X %*% beta)
  sign_eta <- soft_sign(eta, thr)
  
  # Set the initial optimization state
  iter <- 1
  objval <- logpost(y, eta, beta, pL1, pL2)
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
    s  <- sign(eta)*(abs(eta)>thr)
  } else if (type=='PG') {
    w <- c(coeff_PG(eta))
  } else if (type=='NR') { 
    step <- 1.0
    w <- exp(-eta)/(1+exp(-eta))^2
  } else if (type=='BL') {
    w <- rep(0.25, n)
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
  for(iter in 2:maxiter){
    
    # Set the old parameters
    objval_old <- objval
    beta_old <- beta 
    eta_old <- eta
    sign_eta_old <- sign_eta
    
    # Update beta
    beta <- switch(type,
      'NR' = update_enet_coord_NR(y, X, eta, w, beta, pL1, pL2, intercept, niter=1),
      'BL' = update_enet_coord_BL(y, X, eta, w, beta, pL1, pL2, intercept, niter=1),
      'PG' = update_enet_coord_PG(y, X, eta, w, beta, pL1, pL2, intercept, niter=1),
      'PQ' = update_enet_coord_PQ(y, X, eta, w, nu, s, beta, pL1, pL2, intercept, niter=1))
    
    # Update bound
    eta <- as.vector(X %*% beta)
    sign_eta <- soft_sign(eta, thr)
    
    # Check convergence
    objval <- logpost(y, eta, beta, pL1, pL2)
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_beta <- maxreldiff(beta, beta_old)
    diff_rel_eta <- maxreldiff(eta, eta_old)
    min_prod_eta <- min(sign_eta*sign_eta_old)
    
    if (!eta_stable) {
      eta_stable <- (diff_rel_eta<etatol) && (min_prod_eta>=0)
    }
    
    if (type=='PQ'){         
      coeff_pq <- coeff_PQ(eta)
      w  <- c(coeff_pq[,1])
      nu <- c(coeff_pq[,2])
      s  <- (1-phi)*s + phi*sign_eta
      # s  <- ifelse(sign_eta>thr, sign_eta, s)
      # phi <- phi0 * (1 + phi0*iter/10)^-.75
    } else if (type=='PG'){ 
      w <- c(coeff_PG(eta))
    } else if (type=='BL'){
      w <- rep(0.25, n)
    } else if (type=='NR'){
      w <- exp(-eta)/(1+exp(-eta))^2
      # w <- exp(-eta-2*Rmpfr::log1pexp(eta))
      w[(w==0) | is.na(w) | is.infinite(w)] <- 0.25
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

