
#' @title Update beta via PQ bound optimization for ridge logistic regression
#' @keywords internal
update_ridge_PQ <- function(y, X, eta, w, nu, s, pL2, beta, phi, XPXt=NULL,
                            solver=c("chol", "smw", "admm", "sparse"), ctr=NULL){
  # Set the solver
  solver = match.arg(solver)
  # Get the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  XtWy <- crossprod(X, y-0.5-nu*s)
  # Solve the linear system
  # ... ADMM solver
  if (solver=="admm") {
    beta <- admm_genlasso_PQ(y=(y-0.5)/w, X=X, w=w, nu=nu, pL2=pL2, beta0=beta, 
                             lambda=1.0, rho=ctr$rho, precondition=ctr$precondition, 
                             objtol=ctr$objtol, reltol=ctr$reltol, abstol=ctr$abstol, 
                             maxiter=ctr$maxiter, history=FALSE)$beta
  }
  # ... Cholesky solver
  if (solver=="chol") {
    cholQ <- chol(crossprod(X, w*X) + diag(pL2,p,p))
    beta <- chol_solve(cholQ, XtWy)
  }
  # ... Sherman-Morrison-Woodbury solver
  if (solver=="smw") {
    cholQ <- chol(XPXt + diag(1/w,n,n))
    beta <- smw_chol_solve(cholQ, X, pL2, XtWy)
  }
  # ... Sparse Cholesky solver
  if (solver=="sparse") {
    stop("SPARSE-SOLVER: not implemented yet.")
  }
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via PG bound optimization for ridge logistic regression
#' @keywords internal
update_ridge_PG <- function(y, X, eta, w, pL2, XtWy, XPXt=NULL, 
                            solver=c("chol", "smw", "sparse")){
  # Set the solver
  solver = match.arg(solver)
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  if (is.null(XtWy))
    XtWy <- crossprod(X, y-0.5)
  # Solve the linear system
  # ... Cholesky solver
  if (solver=="chol") {
    cholQ <- chol(crossprod(X, w*X) + diag(pL2,p,p))
    beta <- chol_solve(cholQ, XtWy)
  }
  # ... Sherman-Morrison-Woodbury solver
  if (solver=="smw") {
    cholQ <- chol(XPXt + diag(1/w,n,n))
    beta <- smw_chol_solve(cholQ, X, pL2, XtWy)
  }
  # ... Sparse Cholesky solver
  if (solver=="sparse") {
    stop("SPARSE-SOLVER: not implemented yet.")
  }
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via BL bound optimization for ridge logistic regression
#' @keywords internal
update_ridge_BL <- function(y, X, eta, pL2, cholQ,
                            solver=c("chol", "smw", "sparse")){
  # Set the solver
  solver = match.arg(solver)
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  w <- 0.25
  pr <- 1/(1+exp(-eta))
  XtWy <- crossprod(X, y-pr+w*eta)
  # Solve the linear system
  # ... Cholesky solver
  if (solver=="chol")
    beta <- chol_solve(cholQ, XtWy)
  # ... Sherman-Morrison-Woodbury solver
  if (solver=="smw")
    beta <- smw_chol_solve(cholQ, X, pL2, XtWy)
  # ... Sparse Cholesky solver
  if (solver=="sparse")
    stop("SPARSE-SOLVER: not implemented yet.")
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via NR quadratic approximation for ridge logistic regression
#' @keywords internal
update_ridge_NR <- function(y, X, eta, w, pL2, XPXt=NULL, XtWy=NULL,
                            solver=c("chol", "smw", "sparse")){
  # Set the solver
  solver = match.arg(solver)
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the gradient
  pr <- 1/(1+exp(-eta))
  XtWy <- crossprod(X, y-pr+w*eta)
  # Solve the linear system
  # ... Cholesky solver
  if (solver=="chol") {
    cholQ <- chol(crossprod(X, w*X) + diag(pL2,p,p))
    beta <- chol_solve(cholQ, XtWy)
  }
  # ... Sherman-Morrison-Woodbury solver
  if (solver=="smw") {
    cholQ <- chol(XPXt + diag(1/w,n,n))
    beta <- smw_chol_solve(cholQ, X, pL2, XtWy)
  }
  # ... Sparse Cholesky solver
  if (solver=="sparse") {
    stop("SPARSE-SOLVER: not implemented yet.")
  }
  # Return the estimated beta
  return (beta)
}

