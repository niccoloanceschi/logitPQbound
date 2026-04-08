
#' @title Cholesky solvers for dense matrix operations
#' 
#' @param R upper-triangular factor of the Cholesky decomposition
#' @param r vector or matrix on the right-hand side of the linear system
#' @param B cross-product matrix
#'
#' @description 
#' \itemize{
#'   \item \code{chol_solve}: solve a linear system of equations 
#'   \item \code{chol_rnd_sim}: generate a Gaussian random vector 
#'   \item \code{chol_logdet}: log-determinant of a matrix
#'   \item \code{chol_edf}: Hutch++ estimator of the EDF
#' }
#' 
#' @keywords internal
chol_solve <- function(R, r){
  isvec <- is.vector(r)
  ismat <- is.matrix(r)
  r <- forwardsolve(t(R),r)
  r <- backsolve(R,r)
  if(isvec) r <- as.vector(r)
  if(ismat) r <- as.matrix(r)
  return(r)
}

#' @rdname chol_solve
#' @keywords internal
chol_rnd_sim <- function(R, r){
  p <- nrow(R)
  z <- rnorm(p)
  s <- as.vector(forwardsolve(t(R), r))
  x <- as.vector(backsolve(R, z + s))
  return(x)
}

#' @rdname chol_solve
#' @keywords internal
chol_logdet <- function(R){
  logdet <- 2*sum(log(diag(R)))
  return(logdet)
}

#' @rdname chol_solve
#' @keywords internal
chol_edf <- function(R, B, exact=FALSE, rank=10, nsample=100){
  if (exact) {
    edf <- sum(diag(chol_solve(R, B)))
  } else {
    nobs <- nrow(R)
    
    # Randomized range finder
    U <- matrix(rnorm(nobs*(rank+10)), nrow=nobs, ncol=rank+10)
    U <- chol_solve(R, B %*% U)
    Q <- qr.Q(qr(U))[,1:rank]
    V <- crossprod(Q, chol_solve(R, B %*% Q))
    tr_low_rank <- sum(diag(V))
    
    # Residual Hutchinson estimator
    U <- sample(c(-1, +1), size=nobs*nsample, replace=TRUE)
    U <- matrix(U, nrow=nobs, ncol=nsample)
    U <- U - Q %*% crossprod(Q, U)
    V <- B %*% chol_solve(R, U)
    tr_res_hutch <- sum(U * V) / nsample
    
    # Hutch++ estimator
    edf <- tr_low_rank + tr_res_hutch
  }
  
  return(edf)
}

#' @title SMW solvers for dense matrix operations
#' 
#' @param X design matrix
#' @param y response vector
#' @param w weight vector
#' @param pL2 penalty vector
#' @param r vector or matrix on the right-hand side of the linear system
#'
#' @description 
#' \itemize{
#'   \item \code{smw_solve}: solve a linear system of equations 
#'   \item \code{smw_rnd_sim}: generate a Gaussian random vector 
#'   \item \code{smw_logdet}: log-determinant of a matrix
#'   \item \code{smw_edf}: Hutch++ estimator of the EDF
#' }
#' 
#' @keywords internal
smw_solve <- function(X, w, pL2, r){
  isvec <- is.vector(r)
  ismat <- is.matrix(r)
  XP <- t(t(X) / pL2)
  R <- chol(diag(1/w) + tcrossprod(XP, X))
  s <- crossprod(XP, chol_solve(R, XP %*% r))
  if(isvec) s <- as.vector(s)
  if(ismat) s <- as.matrix(s)
  return(r/pL2 - s)
}

#' @rdname smw_solve
#' @keywords internal
smw_rnd_sim <- function(X, y, pL2){
  n <- nrow(X)
  p <- ncol(X)
  u <- rnorm(p, mean=0, sd=sqrt(1/pL2))
  d <- rnorm(n, mean=0, sd=1)
  v <- as.vector(X %*% u + d)
  w <- as.vector(solve(X %*% (t(X)/pL2) + diag(n), y-v))
  z <- u + as.vector(crossprod(X, w)/pL2)
  return(z)
}

#' @rdname smw_solve
#' @keywords internal
smw_logdet <- function(X, w, pL2){
  R <- chol(diag(1/w) + tcrossprod(t(t(X)/sqrt(pL2))))
  logdet <- sum(log(pL2)) + sum(log(w)) + 2*sum(log(diag(R)))
  return(logdet)
}

#' @rdname smw_solve
#' @keywords internal
smw_edf <- function(X, w, pL2){
  n <- nrow(X)
  p <- ncol(X)
  XPXt <- X %*% (t(X) / pL2)
  G <- XPXt %*% (diag(n) - solve(diag(1/w) + XPXt, XPXt)) %*% diag(w)
  edf <- sum(diag(G))
  return(edf)
}

#' @title SMW-Cholesky solvers for dense matrix operations
#' 
#' @param R upper-triangular factor of the Cholesky decomposition
#' @param X design matrix
#' @param B cross-product matrix
#' @param w weight vector
#' @param pL2 penalty vector
#' @param r vector or matrix on the right-hand side of the linear system
#'
#' @description 
#' \itemize{
#'   \item \code{smw_chol_solve}: solve a linear system of equations 
#'   \item \code{smw_chol_rnd_sim}: generate a Gaussian random vector 
#'   \item \code{smw_chol_logdet}: log-determinant of a matrix
#'   \item \code{smw_chol_edf}: Hutch++ estimator of the EDF
#' }
#' 
#' @keywords internal
smw_chol_solve <- function(R, X, pL2, r){
  isvec <- is.vector(r)
  ismat <- is.matrix(r)
  s <- (1/pL2) * crossprod(X, chol_solve(R, X %*% (r/pL2)))
  if(isvec) s <- as.vector(s)
  if(ismat) s <- as.matrix(s)
  return(r/pL2 - s)
}

#' @rdname smw_chol_solve
#' @keywords internal
smw_chol_logdet <- function(R, w, pL2){
  logdet <- sum(log(pL2)) + sum(log(w)) + 2*sum(log(diag(R)))
  return(logdet)
}

#' @rdname smw_chol_solve
#' @keywords internal
smw_chol_edf <- function(R, B, w, pL2){
  n <- nrow(X)
  p <- ncol(X)
  G <- B %*% (diag(n) - chol_solve(R, B)) %*% diag(w)
  edf <- sum(diag(G))
  return(edf)
}


#' @title Cholesky solvers for sparse matrix operations
#' 
#' @param R upper-triangular factor of the Cholesky decomposition
#' @param r vector or matrix on the right-hand side of the linear system
#' @param B cross-product matrix
#' @param exact if \code{TRUE}, compute the exact edf
#' @param rank number of eigen-vectors to use
#' @param nsample number of random sample
#' @param seed random number generation seed
#' 
#' @description 
#' \itemize{
#'   \item \code{sp_chol_solve}: solve a linear system of equations 
#'   \item \code{sp_chol_rnd_sim}: generate a Gaussian random vector 
#'   \item \code{sp_chol_logdet}: log-determinant of a matrix
#'   \item \code{sp_chol_edf}: Hutch++ estimator of the EDF
#' }
#' 
#' @keywords internal
sp_chol_solve <- function(R, r){
  isvec <- is.vector(r)
  ismat <- is.matrix(r) | is(r, "Matrix")
  r <- Matrix::solve(R,r)
  if(isvec) r <- as.vector(r)
  if(ismat) r <- as.matrix(r)
  return(r)
}

#' @rdname sp_chol_solve
#' @keywords internal
sp_chol_rnd_sim <- function(R, r){
  p <- nrow(R)
  z <- rnorm(p)
  s <- as.vector(Matrix::solve(R, r, system="L")) # forwardsolve
  x <- as.vector(Matrix::solve(R, z+s, system="Lt")) # backsolve
  return(x)
}

#' @rdname sp_chol_solve
#' @keywords internal
sp_chol_logdet <- function(R){
  logdet <- 2*sum(log(Matrix::diag(Matrix::expand(R)$L)))
  return(logdet)
}

#' @rdname sp_chol_solve
#' @keywords internal
sp_chol_edf <- function(R, B, exact=FALSE, rank=10, nsample=100, seed=1234){
  if (exact) {
    # See also the package sparseinv
    edf <- sum(Matrix::diag(Matrix::solve(R, B)))
  } else {
    set.seed(seed)
    
    # Number of observations
    nobs <- nrow(R)
    
    # Randomized range finder
    U <- matrix(rnorm(nobs*(rank+10)), nrow=nobs, ncol=rank+10)
    U <- sp_chol_solve(R, B %*% U)
    Q <- qr.Q(qr(U))[,1:rank]
    V <- crossprod(Q, sp_chol_solve(R, B %*% Q))
    tr_low_rank <- sum(diag(V))
    
    # Residual Hutchinson estimator
    U <- sample(c(-1, +1), size=nobs*nsample, replace=TRUE)
    U <- matrix(U, nrow=nobs, ncol=nsample)
    U <- as.matrix(U - Q %*% crossprod(Q, U))
    V <- as.matrix(B %*% sp_chol_solve(R, U))
    tr_res_hutch <- sum(U * V) / nsample
    
    # Hutch++ estimator
    edf <- tr_low_rank + tr_res_hutch
  }
  
  return(edf)
}

#' @title Maximum absolute relative difference operator
#' 
#' @param xnew new scalar, vector or matrix value
#' @param xold old scalar, vector or matrix value
#' 
#' @keywords internal
maxreldiff <- function(xnew, xold) {
  return(max(abs(xnew-xold) / (abs(xold)+1e-8)))
}

#' @title Element-wise soft-sign operator
#' 
#' @param x vector or matrix to be thresholded
#' @param thr thresholding parameter
#' 
#' @keywords internal
soft_sign <- function(x, thr) {
  return(sign(x) * (abs(x) > thr))
}

#' @title Element-wise soft-thresholding operator
#' 
#' @param x vector or matrix to be thresholded
#' @param lambda thresholding parameter
#' 
#' @keywords internal
soft_threshold <- function(x, lambda) {
  return(sign(x) * pmax(0, abs(x)-lambda))
}

