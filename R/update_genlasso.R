
#' @title Initialize the regression parameters for generalized-lasso logistic regression
#' @keywords internal
init_genlasso_beta = function(y, X, D, pL1, pL2, spthr=0.9){
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the parameter vector
  beta <- numeric(p)
  
  # Check for sparsity
  if (sum(abs(D)<1e-8) > spthr*m*p) {
    D[abs(D) < 1e-8] <- 0
    D = Matrix::Matrix(D, sparse=TRUE)
  }
  if (sum(abs(X) < 1e-8) > spthr*n*p) {
    X[abs(X)<1e-8] <- 0
    X = Matrix::Matrix(X, sparse=TRUE)
  }
  
  # Compute the penalty matrix
  P <- Matrix::crossprod(D, pL1 * D) + diag(pL2)
  
  if (p>n) {
    # If p>n, we use the SMW identity + the sparsity of P
    cholP <- Matrix::Cholesky(P, perm=TRUE)
    PXty <- sp_chol_solve(cholP, Matrix::crossprod(X, y-0.5))
    XP <- t(sp_chol_solve(cholP, t(as.matrix(X))))
    Q <- Matrix::tcrossprod(X, XP) + diag(4,n,n)
    cholQ <- base::chol(as.matrix(Q))
    beta[] <- PXty - as.vector(Matrix::crossprod(XP, chol_solve(cholQ, X %*% PXty)))
  } else {
    # If p<n, we use direct inversion via Cholesky factorization
    Q <- 0.25 * Matrix::crossprod(X) + P
    if (sum(abs(Q)<1e-8) > spthr*p*p) {
      Q[abs(Q)<1e-8] <- 0
      Q <- Matrix::Matrix(Q, sparse=TRUE)
      cholQ <- Matrix::Cholesky(Q, perm=TRUE)
      beta[] <- sp_chol_solve(cholQ, crossprod(X, y-0.5))
    } else {
      cholQ <- base::chol(as.matrix(Q))
      beta[] <- chol_solve(cholQ, crossprod(X, y-0.5))
    }
  }
  
  return(beta)
}

#' @title Update beta via PQ bound optimization for generalized-lasso logistic regression
#' @keywords internal
update_genlasso_PQ = function(beta, z, u, s, eta, y, X, w, nu, D, lambda, alpha, eps, phi, approx, ctr){
  # stable <- FALSE
  # z <- NULL
  if (!approx) {
    # If sign(eta) did not stabilize yet, run the augmented ADMM algorithm
    r <- (y - 0.5) / w
    beta <- solve_admm_genenet_PQ(
      X=X, y=r, w=w, nu=nu, D=D, beta0=beta, z0=z, u0=u, lambda=lambda, alpha=alpha, 
      eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
      spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
      objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
      maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  } else {
    # If sign(eta) is stable, run the reduced ADMM algorithm
    n <- length(y)
    # r <- (y - 0.5 - nu * sign(eta)) / w
    thr <- 1e-2
    s <- (1-phi)*s + phi*(sign(eta)*(abs(eta)>thr))
    r <- (y - 0.5 - nu*s) / w
    # if (length(z)>nrow(D)) {
    #   z0 <- z[1:n]
    #   u0 <- u[1:n]
    #   z <- z[-(1:n)]
    #   u <- u[-(1:n)]
    # } 
    beta <- solve_admm_genenet_WLS(
      X=X, y=r, w=w, D=D, beta0=beta, z0=z, u0=u, lambda=lambda, alpha=alpha, 
      eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
      spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
      objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
      maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
    # 
    # beta$z <- c(z0, beta$z)
    # beta$u <- c(u0, beta$u)
    # Rename the primal and dual residuals
    beta$r_pri  <- beta$r
    beta$r_dual <- beta$s
    # Store the pseudo-data and smoothed signs
    beta$r <- r
    beta$s <- s
  }
  return(beta)
}

#' @title Update beta via PG bound optimization for generalized-lasso logistic regression
#' @keywords internal
update_genlasso_PG = function(beta, z, u, y, X, w, D, lambda, alpha, eps, ctr){
  z <- NULL
  beta <- solve_admm_genenet_WLS(
    X=X, y=(y-0.5)/w, w=w, D=D, beta0=beta, z0=z, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}

#' @title Update beta via BL bound optimization for generalized-lasso logistic regression
#' @keywords internal
update_genlasso_BL = function(beta, z, u, eta, prob, y, X, D, lambda, alpha, eps, ctr){
  # WARNING: here cholQ is recomputed inside any instance of the ADMM alg, this must be fixed!
  z <- NULL
  w <- rep(0.25, length(y))
  r <- (y - prob + w * eta) / w
  beta <- solve_admm_genenet_WLS(
    X=X, y=r, w=w, D=D, beta0=beta, z0=z, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}

#' @title Update beta via NR quadratic approximation for generalized-lasso logistic regression
#' @keywords internal
update_genlasso_NR = function(beta, z, u, eta, prob, y, X, w, D, lambda, alpha, eps, ctr){
  z <- NULL
  w <- ifelse(w<1e-8, 1e-8, w) 
  r <- (y - prob + w * eta) / w
  beta <- solve_admm_genenet_WLS(
    X=X, y=r, w=w, D=D, beta0=beta, z0=z, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}



#' @title Initialize the regression parameters for spatial-lasso logistic regression
#' @keywords internal
init_splasso_beta = function(y, X, D, pL1, pL2, spthr=0.9){
  
  # Set the parameter dimension
  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(D)
  
  # Compute the sufficient statistics
  Q <- Matrix::crossprod(rbind(0.5*X, sqrt(pL1)*D, diag(sqrt(pL2),p,p)))
  R <- Matrix::Cholesky(Q, perm=TRUE)
  r <- as.vector(Matrix::crossprod(X, y-0.5))
  
  # Initialize beta
  beta <- sp_chol_solve(R, r)
  
  # Output
  return(beta)
}


#' @title Update beta via PQ bound optimization for spatial-lasso logistic regression
#' @keywords internal
update_splasso_PQ = function(beta, z, u, s, eta, y, X, w, nu, D, lambda, alpha, eps, phi, approx, ctr){
  if (!approx) {
    # EXACT optimization: run the augmented ADMM algorithm
    beta <- solve_admm_spenet_PQ(
      X=X, y=(y-0.5)/w, w=w, nu=nu, D=D, beta0=beta, z0=NULL, u0=u, 
      lambda=lambda, alpha=alpha, eps=eps, rho=ctr$rho, 
      tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
      spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
      objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
      maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  } else {
    # APPROXIMATE optimization: run the reduced ADMM algorithm with active-set correction
    thr <- 1e-2
    s <- (1-phi)*s + phi*(sign(eta)*(abs(eta)>thr))
    r <- (y - 0.5 - nu*s) / w
    beta <- solve_admm_spenet_WLS(
      X=X, y=r, w=w, D=D, beta0=beta, z0=z, u0=u, 
      lambda=lambda, alpha=alpha, eps=eps, rho=ctr$rho, 
      tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
      spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
      objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
      maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
    # Rename the primal and dual residuals
    beta$r_pri  <- beta$r
    beta$r_dual <- beta$s
    # Store the pseudo-data and smoothed signs
    beta$r <- r
    beta$s <- s
  }
  return(beta)
}

#' @title Update beta via PG bound optimization for spatial-lasso logistic regression
#' @keywords internal
update_splasso_PG = function(beta, z, u, y, X, w, D, lambda, alpha, eps, ctr){
  beta <- solve_admm_spenet_WLS(
    X=X, y=(y-0.5)/w, w=w, D=D, beta0=beta, z0=NULL, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}

#' @title Update beta via BL bound optimization for spatial-lasso logistic regression
#' @keywords internal
update_splasso_BL = function(beta, z, u, eta, prob, y, X, D, lambda, alpha, eps, ctr){
  # WARNING: here cholQ is recomputed inside any instance of the ADMM alg, this must be fixed!
  w <- rep(0.25, length(y))
  r <- (y - prob + w * eta) / w
  beta <- solve_admm_spenet_WLS(
    X=X, y=r, w=w, D=D, beta0=beta, z0=NULL, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}

#' @title Update beta via NR quadratic approximation for spatial-lasso logistic regression
#' @keywords internal
update_splasso_NR = function(beta, z, u, eta, prob, y, X, w, D, lambda, alpha, eps, ctr){
  w <- ifelse(w<1e-8, 1e-8, w) 
  r <- (y - prob + w * eta) / w
  beta <- solve_admm_spenet_WLS(
    X=X, y=r, w=w, D=D, beta0=beta, z0=NULL, u0=u, lambda=lambda, alpha=alpha, 
    eps=eps, rho=ctr$rho, tau=ctr$tau, gamma=ctr$gamma, intercept=ctr$intercept, 
    spthr=ctr$spthr, smw=ctr$smw, precondition=ctr$precondition, 
    objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
    maxiter=ctr$maxiter, verbose=ctr$verbose, freq=ctr$freq)
  
  return(beta)
}

