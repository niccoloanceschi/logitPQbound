
#' ADMM optimizer for the PQ generalized lasso problem
#' 
#' @param X design matrix
#' @param y response vector
#' @param w weight vector
#' @param nu penalty weight vector
#' @param pL2 ridge penalty vector
#' @param beta0 initial values
#' @param lambda generalized lasso penalty parameter
#' @param rho augmented Lagrangian penalty parameter
#' @param objtol objective tolerance parameter
#' @param reltol relative tolerance parameter
#' @param abstol absolute tolerance parameter
#' @param maxiter maximum number of iterations
#' 
#' @return A list containing the estimated coefficients and the optimization history.
#' 
#' @keywords internal
solve_admm_ridge_PQ <- function(X, y, w, nu, pL2, beta0=NULL, 
                                rho=1.0, precondition=FALSE, 
                                objtol=1e-5, reltol=1e-2, abstol=1e-4, 
                                maxiter=1000, history=FALSE) {
  
  # Set the initial values
  if (is.null(beta0)) {
    p <- ncol(X)
    beta0 <- rnorm(p, mean=0, sd=sqrt(1/p))
    beta0 <- beta0 / sqrt(sum(beta0^2))
  }
  
  # Call the Rcpp implementation
  out <- admm_genlasso_PQ(y=y, X=X, w=w, nu=nu, pL2=pL2, beta0=beta0, 
                          lambda=1.0, rho=rho, precondition=precondition, 
                          objtol=objtol, reltol=reltol, abstol=abstol, 
                          maxiter=maxiter, history=history)
  
  # transform the trace into a data-frame
  if (history) {
    output$trace <- as.data.frame(output$trace)
    colnames(output$trace) = c("iter", "objval", "norm_b", "norm_z", "norm_u", 
                               "norm_r", "norm_s", "eps_pri", "eps_dual")
  }
  
  # Return the estimated object
  return(out)
}

#' Compute the generalized elastic-net solution using the ADMM algorithm
#' 
#' @param X design matrix
#' @param y response vector
#' @param w L2 weight vector 
#' @param nu L1 weight vector 
#' @param D penalty matrix
#' @param beta0 optional initial value for beta
#' @param z0 optional initial value for the primal variables
#' @param u0 optional initial value for the dual variables
#' @param lambda penalty parameter
#' @param alpha elastic-net mixing parameter 
#' @param eps intercept specific penalty parameter
#' @param rho augmented Lagrangian penalty parameter
#' @param tau dual over-relaxation parameter
#' @param gamma proximal regularization parameter
#' @param intercept if `TRUE`, do not penalize per first column of `X` 
#' @param spthr proportion of zero entries to transform `X` and `D` into sparse matrices
#' @param smw if `TRUE`, allows for the SMW solver in the `p>n` regime 
#' @param precondition if `TRUE`, allows for diagonal preconditioning
#' @param objtol objective tolerance stopping criterion
#' @param reltol relative tolerance stopping criterion
#' @param abstol absolute tolerance stopping criterion
#' @param maxiter maximum number of iterations
#' @param verbose if `TRUE`, print the current optimization status
#' @param freq frequency of printing of the optimization status
#' 
#' @return 
#' \itemize{
#'   \item \code{beta}: estimated primal coefficients 
#'   \item \code{z}: estimated primal variable
#'   \item \code{u}: estimated dual variables
#'   \item \code{r}: primal residuals
#'   \item \code{s}: dual residuals
#'   \item \code{niter}: number of iterations to reach convergence
#'   \item \code{trace}: data-frame collecting the optimization history
#' }
#' 
#' @keywords internal
solve_admm_genenet_PQ <- function(X, y, w, nu, D=diag(ncol(X)), d=numeric(nrow(D)),
                                  beta0=NULL, z0=NULL, u0=NULL, lambda=1.0, alpha=0.99, 
                                  eps=1e-8, rho=1.0, tau=0.0, gamma=0.0, intercept=TRUE, 
                                  spthr=0.9, smw=FALSE, precondition=FALSE,
                                  objtol=1e-5, reltol=1e-2, abstol=1e-4, 
                                  maxiter=1000, verbose=FALSE, freq=10) {
  
  # Model dimensions
  m <- nrow(D)
  n <- nrow(X)
  p <- ncol(X)
  
  sqrtm <- sqrt(m)
  sqrtn <- sqrt(n)
  sqrtp <- sqrt(p)
  sqrtnm <- sqrt(n+m)
  
  # Check for sparsity
  is_sp_D = FALSE
  is_sp_X = FALSE
  is_sp_P = FALSE
  is_sp_Q = FALSE
  
  if (sum(abs(D)<1e-8) > spthr*m*p) {
    D[abs(D) < 1e-8] <- 0
    D = Matrix::Matrix(D, sparse=TRUE)
    is_sp_D = TRUE
  }
  if (sum(abs(X) < 1e-8) > spthr*n*p) {
    X[abs(X)<1e-8] <- 0
    X = Matrix::Matrix(X, sparse=TRUE)
    is_sp_X = TRUE
  }
    
  # Set the penalty vectors
  if ((alpha>1) || (alpha<=0)) {
    stop("ADMM-PQ-ENET: alpha must be real value in the interval (0,1].")
  }
  
  pL1 <- rep(lambda*alpha, m)
  pL2 <- rep(lambda*(1-alpha), p)
  if (intercept) {
    pL2[1] <- eps*lambda*(1-alpha)
  }
  
  # Preconditioning
  pwts <- rep(1, times=n+m)
  if (precondition) {
    pwts[ (1:n)] <- 1 / sqrt(Matrix::rowSums(X^2) * nu^2 + 1e-6)
    pwts[-(1:n)] <- 1 / sqrt(Matrix::rowSums(D^2) + 1e-6)
    
    pwts[ (1:n)] <- 1 / (Matrix::rowSums(abs(X)) * abs(nu) + 1e-6)
    pwts[-(1:n)] <- 1 / (Matrix::rowSums(abs(D)) + 1e-6)
  }
    
  # Initialization
  # ... Primal variables: regression coefficients
  beta <- if (!is.null(beta0)) beta0 else rnorm(p, sd=0.1)
  # ... Primal variables: linear predictor
  eta <- as.vector(X %*% beta)
  NXb <- as.vector(nu * eta)
  Db <- as.vector(D %*% beta)
  Hb <- c(NXb, Db - d)
  # Hb <- c(NXb, lambda*alpha*(Db - d))
  # ... Primal variables: augmented variables
  z <- if (!is.null(z0)) z0 else Hb
  # ... Dual variables: Lagrange multipliers
  u <- if (!is.null(u0)) u0 else numeric(n+m) 
  # ... Primal variables: previous iteration
  z_old <- z
  # ... Primal and dual residuals
  r <- numeric(n+m)
  s <- numeric(n+m)
  
  # Pre-compute static quantities
  # ... Working weights, weighted responses, multiplicative factors
  wts <- w + rho*nu^2*pwts[1:n]
  d0 <- c(rep(0, n), d)
  wy_0 <- c(w * y, rep(0, m))
  # rho_nu_0 <- rho * c(nu, rep(lambda*alpha, m))
  rho_nu_0 <- rho * pwts * c(nu, rep(1, m))
  pen_0 <- c(rep(1/rho, n), pL1/rho) / pwts
  # ... System matrix factorization
  # P <- rho*(lambda*alpha)^2 * Matrix::crossprod(D) + diag(pL2)
  P <- rho * Matrix::crossprod(D, pwts[-(1:n)]*D) + diag(pL2 + gamma)
  if (p>n && smw) {
    # If p>n, we use the SMW identity + the sparsity of P
    cholP <- Matrix::Cholesky(P, pivot=TRUE)
    XP <- t(sp_chol_solve(cholP, t(X)))
    Q <- Matrix::tcrossprod(X, XP) + diag(1/wts,n,n)
    cholQ <- base::chol(as.matrix(Q))
    is_sp_Q <- FALSE
  } else {
    # If p<n, we use sparse pivoted Cholesky factorization
    Q <- Matrix::crossprod(X, wts*X) + P
    if (sum(abs(Q)<1e-8) > spthr*p*p) {
      Q[abs(Q)<1e-8] <- 0
      Q <- Matrix::Matrix(Q, sparse=TRUE)
      cholQ <- Matrix::Cholesky(Q, pivot=TRUE)
      is_sp_Q <- TRUE
    } else {
      cholQ <- base::chol(as.matrix(Q))
      is_sp_Q <- FALSE
    }
  }
  
  # Objective function 
  loss_PQ <- 0.5 * sum(w*(y - eta)^2) + sum(abs(NXb))
  pen_ridge <- 0.5 * sum(pL2*beta^2)
  pen_lasso <- sum(pL1*abs(Db - d))
  objval <- loss_PQ + pen_ridge + pen_lasso
  objval_old <- objval
  objval_min <- objval
  
  # Set the optimization history
  trace <- as.data.frame(matrix(NaN, nrow=maxiter, ncol=10))
  colnames(trace) <- c("iter", "objval", "diffobj", "norm_b", "norm_z", "norm_u", 
                       "norm_r", "norm_s", "eps_pri", "eps_dual")
  
  # Store the initial values
  trace$iter[1] <- 1
  trace$objval[1] <- objval
  trace$diffobj[1] <- NaN
  trace$norm_b[1] <- sqrt(sum(beta^2)) / sqrtp
  trace$norm_z[1] <- sqrt(sum(z^2)) / sqrtnm
  trace$norm_u[1] <- sqrt(sum(u^2)) / sqrtnm
  trace$norm_r[1] <- NaN
  trace$norm_s[1] <- NaN
  trace$eps_pri[1] <- NaN
  trace$eps_dual[1] <- NaN
  
  # Initial output
  iter = 1
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "\n")
  }
  
  # ADMM optimization loop
  for (iter in 2:maxiter) {
    
    # Primal update: beta
    e <- as.vector(Matrix::crossprod(rbind(X, D), wy_0 + rho_nu_0 * (d0+z-u))) + gamma*beta
    if (p>n && smw) {
      f <- sp_chol_solve(cholP, e)
      beta[] <- f - as.vector(Matrix::crossprod(XP, chol_solve(cholQ, X %*% f)))
    } else {
      if (is_sp_Q) {
        beta[] <- sp_chol_solve(cholQ, e)
      } else {
        beta[] <- chol_solve(cholQ, e)
      }
    }
    
    # Primal update: eta, NXb, Db
    eta[] <- as.vector(X %*% beta)
    NXb[] <- as.vector(nu * eta)
    Db[] <- as.vector(D %*% beta)
    Hb[] <- c(NXb, Db - d)
    # Hb[] <- c(NXb, lambda*alpha*(Db - d))
    
    # I dual update: u
    u[] <- u + tau * (Hb - z)
    
    # Primal update: z
    z_old[] <- z
    z[] <- soft_threshold(c(Hb + u), pen_0)

    # II dual update: u
    u[] <- u + (1-tau) * (Hb - z)

    # Primal and dual residuals
    r[] <- Hb - z
    s[] <- - rho * pwts * (z - z_old)

    # Convergence diagnostics
    norm_b <- sqrt(sum(beta^2))
    norm_z <- sqrt(sum(z^2))
    norm_u <- sqrt(sum(u^2))
    norm_r <- sqrt(sum(r^2))
    norm_s <- sqrt(sum(s^2))
    eps_pri <- max(sqrtp, sqrtnm) * abstol + reltol * max(norm_b, norm_z)
    eps_dual <- sqrtnm * abstol + reltol* rho * norm_u
    
    # Objective function
    loss_PQ <- 0.5 * sum(w*(y - eta)^2) + sum(abs(NXb))
    pen_ridge <- 0.5 * sum(pL2*beta^2)
    pen_lasso <- sum(pL1*abs(Db - d))
    objval <- loss_PQ + pen_ridge + pen_lasso
    diffobj <- abs(objval-objval_old) / (abs(objval)+1e-8)
    objval_min <- min(objval_min, objval_old)
    objval_old <- objval
    
    # Optimization state allocation
    trace$iter[iter] <- iter
    trace$objval[iter] <- objval
    trace$diffobj[iter] <- diffobj
    trace$norm_b[iter] <- norm_b / sqrtp
    trace$norm_z[iter] <- norm_z / sqrtm
    trace$norm_u[iter] <- norm_u / sqrtm
    trace$norm_r[iter] <- norm_r / sqrtm
    trace$norm_s[iter] <- norm_s / sqrtm
    trace$eps_pri[iter] <- eps_pri
    trace$eps_dual[iter] <- eps_dual
    
    # Intermediate output
    if (verbose && iter%%freq==0) {
      cat("  iter:", gettextf("%3d", iter),
          "  objval:", gettextf("%.4f", objval/n),
          "  diffobj:", gettextf("%.5f", diffobj),
          "  norm_r:", gettextf("%.4f", norm_r),
          "  norm_s:", gettextf("%.4f", norm_s),
          "\n")
    }
    
    # Convergence check
    chech_obj_sign <- (objval < objval_min)
    check_obj_diff <- (diffobj < objtol)
    check_r_norm <- (norm_r < eps_pri)
    check_s_norm <- (norm_s < eps_dual)
    if (chech_obj_sign && check_obj_diff && check_r_norm && check_s_norm) {
      break
    }
  }
  
  # Remove unused iterations
  if (iter < maxiter){
    trace <- trace[1:iter,]
  }
  
  # Return the results
  return(list(beta=beta, z=z, u=u, r=r, s=s, niter=iter, trace=trace))
}

#' @rdname solve_admm_genenet_PQ
#' @keywords internal
solve_admm_genenet_WLS <- function(X, y, w, D=diag(ncol(X)), d=numeric(nrow(D)),
                                   beta0=NULL, z0=NULL, u0=NULL, lambda=1.0, alpha=0.99, 
                                   eps=1e-8, rho=1.0, tau=0.0, gamma=0.0, intercept=TRUE, 
                                   spthr=0.9, smw=FALSE, precondition=FALSE,
                                   objtol=1e-5, reltol=1e-2, abstol=1e-4, 
                                   maxiter=1000, verbose=FALSE, freq=10) {
  
  # Model dimensions
  m <- nrow(D)
  n <- nrow(X)
  p <- ncol(X)
  
  sqrtm <- sqrt(m)
  sqrtn <- sqrt(n)
  sqrtp <- sqrt(p)
  
  # Check for sparsity
  is_sp_D = FALSE
  is_sp_X = FALSE
  is_sp_P = FALSE
  is_sp_Q = FALSE
  
  if (sum(abs(D)<1e-8) > spthr*m*p) {
    D[abs(D) < 1e-8] <- 0
    D = Matrix::Matrix(D, sparse=TRUE)
    is_sp_D = TRUE
  }
  if (sum(abs(X) < 1e-8) > spthr*n*p) {
    X[abs(X)<1e-8] <- 0
    X = Matrix::Matrix(X, sparse=TRUE)
    is_sp_X = TRUE
  }
  
  # Set the penalty vectors
  if ((alpha>1) || (alpha<=0)) {
    stop("ADMM-WLS-ENET: alpha must be real value in the interval (0,1].")
  }
  
  pL1 <- rep(lambda*alpha, m)
  pL2 <- rep(lambda*(1-alpha), p)
  if (intercept) {
    pL2[1] <- eps*lambda*(1-alpha)
  }
  
  # Preconditioning
  if (precondition) {
    pwts <- 1 / sqrt(Matrix::rowSums(D^2)+1e-6)
    pwts <- 1 / (Matrix::rowSums(abs(D)) + 1e-6)
  } else {
    pwts <- rep(1, times=m)
  }
  
  # Initialization
  # ... Primal variables: regression coefficients
  beta <- if (!is.null(beta0)) beta0 else rnorm(p, sd=0.1)
  # ... Primal variables: linear predictor
  eta <- as.vector(X %*% beta)
  Db <- as.vector(D %*% beta)
  # ... Primal variables: augmented variables
  z <- if (!is.null(z0)) z0 else (Db - d)
  # ... Dual variables: Lagrange multipliers
  u <- if (!is.null(u0)) u0 else numeric(m)
  # ... Primal variables: previous iteration
  z_old <- z
  # ... Primal and dual residuals
  r <- numeric(m)
  s <- numeric(m)
  
  # Pre-compute static quantities
  # ... Working weights, weighted responses, multiplicative factors
  XWy <- as.vector(Matrix::crossprod(X, w * y))
  # ... System matrix factorization
  # P <- rho*(lambda*alpha)^2 * Matrix::crossprod(D) + diag(pL2)
  P <- rho * Matrix::crossprod(D, pwts*D) + diag(pL2 + gamma)
  if (p>n && smw) {
    # If p>n, we use the SMW identity + the sparsity of P
    cholP <- Matrix::Cholesky(P, perm=TRUE)
    XP <- t(sp_chol_solve(cholP, t(as.matrix(X))))
    Q <- Matrix::tcrossprod(X, XP) + diag(1/w,n,n)
    cholQ <- base::chol(as.matrix(Q))
    is_sp_Q <- FALSE
  } else {
    # If p<n, we use direct inversion via Cholesky factorization
    Q <- Matrix::crossprod(X, w*X) + P
    if (sum(abs(Q)<1e-8) > spthr*p*p) {
      Q[abs(Q)<1e-8] <- 0
      Q <- Matrix::Matrix(Q, sparse=TRUE)
      cholQ <- Matrix::Cholesky(Q, perm=TRUE)
      is_sp_Q <- TRUE
    } else {
      cholQ <- base::chol(as.matrix(Q))
      is_sp_Q <- FALSE
    }
  }
  
  # Objective function 
  loss_WLS <- 0.5 * sum(w*(y - eta)^2)
  pen_ridge <- 0.5 * sum(pL2*beta^2)
  pen_lasso <- sum(pL1*abs(Db - d))
  objval <- loss_WLS + pen_ridge + pen_lasso
  objval_old <- objval
  objval_min <- objval
  
  # Set the optimization history
  trace <- as.data.frame(matrix(NaN, nrow=maxiter, ncol=10))
  colnames(trace) <- c("iter", "objval", "diffobj", "norm_b", 
                       "norm_z", "norm_u", "norm_r", "norm_s", 
                       "eps_pri", "eps_dual")
  
  # Store the initial values
  trace$iter[1] <- 1
  trace$objval[1] <- objval
  trace$diffobj[1] <- NaN
  trace$norm_b[1] <- sqrt(sum(beta^2)) / sqrtp
  trace$norm_z[1] <- sqrt(sum(z^2)) / sqrtm
  trace$norm_u[1] <- sqrt(sum(u^2)) / sqrtm
  trace$norm_r[1] <- NaN
  trace$norm_s[1] <- NaN
  trace$eps_pri[1] <- NaN
  trace$eps_dual[1] <- NaN
  
  # Initial output
  iter = 1
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "\n")
  }
  
  # ADMM optimization loop
  for (iter in 2:maxiter) {
    
    # Primal update: beta
    # e <- XWy + (rho*lambda*alpha) * as.vector(Matrix::crossprod(D, d+z-u))
    e <- XWy + rho * as.vector(Matrix::crossprod(D, pwts*(d+z-u))) + gamma*beta
    if (p>n && smw) {
      f <- sp_chol_solve(cholP, e)
      beta[] <- f - as.vector(Matrix::crossprod(XP, chol_solve(cholQ, X %*% f)))
    } else {
      if (is_sp_Q) {
        beta[] <- sp_chol_solve(cholQ, e)
      } else {
        beta[] <- chol_solve(cholQ, e)
      }
    }
    
    # Primal update: eta, NXb, Db
    eta[] <- as.vector(X %*% beta)
    Db[] <- as.vector(D %*% beta)
    
    # I dual update: u 
    u[] <- u + tau * (Db - d - z)
    
    # Primal update: z
    z_old[] <- z
    z[] <- soft_threshold(c(Db - d + u), (lambda*alpha/rho)/pwts)
    
    # II dual update: u
    u[] <- u + (1-tau) * (Db - d - z)
    
    # Primal and dual residuals
    r[] <- Db - d - z
    s[] <- - rho * pwts * (z - z_old)
    
    # Convergence diagnostics
    norm_b <- sqrt(sum(beta^2))
    norm_z <- sqrt(sum(z^2))
    norm_u <- sqrt(sum(u^2))
    norm_r <- sqrt(sum(r^2))
    norm_s <- sqrt(sum(s^2))
    eps_pri <- abstol * max(sqrtp, sqrtm) + reltol * max(norm_b, norm_z)
    eps_dual <- abstol * sqrtm + reltol * rho * norm_u
    
    # Objective function
    loss_WLS <- 0.5 * sum(w*(y - eta)^2)
    pen_ridge <- 0.5 * sum(pL2*beta^2)
    pen_lasso <- sum(pL1*abs(Db - d))
    objval <- loss_WLS + pen_ridge + pen_lasso
    diffobj <- abs(objval-objval_old) / (abs(objval)+1e-8)
    objval_min <- min(objval_min, objval_old)
    objval_old <- objval
    
    # Optimization state allocation
    trace$iter[iter] <- iter
    trace$objval[iter] <- objval
    trace$diffobj[iter] <- diffobj
    trace$norm_b[iter] <- norm_b / sqrtp
    trace$norm_z[iter] <- norm_z / sqrtm
    trace$norm_u[iter] <- norm_u / sqrtm
    trace$norm_r[iter] <- norm_r / sqrtm
    trace$norm_s[iter] <- norm_s / sqrtm
    trace$eps_pri[iter] <- eps_pri
    trace$eps_dual[iter] <- eps_dual
    
    # Intermediate output
    if (verbose && iter%%freq==0) {
      cat("  iter:", gettextf("%3d", iter),
          "  objval:", gettextf("%.4f", objval/n),
          "  diffobj:", gettextf("%.5f", diffobj),
          "  norm_r:", gettextf("%.4f", norm_r),
          "  norm_s:", gettextf("%.4f", norm_s),
          "\n")
    }
    
    # Convergence check
    chech_obj_sign <- (objval < objval_min)
    check_obj_diff <- (diffobj < objtol)
    check_r_norm <- (norm_r < eps_pri)
    check_s_norm <- (norm_s < eps_dual)
    if (chech_obj_sign && check_obj_diff && check_r_norm && check_s_norm) {
      break
    }
  }
  
  # Intermediate output
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "  diffobj:", gettextf("%.5f", diffobj),
        "  norm_r:", gettextf("%.4f", norm_r),
        "  norm_s:", gettextf("%.4f", norm_s),
        "\n")
  }
  
  # Remove unused iterations
  if (iter < maxiter){
    trace <- trace[1:iter,]
  }
  
  # Return the results
  return(list(beta=beta, z=z, u=u, r=r, s=s, niter=iter, trace=trace))
}


#' Compute the generalized elastic-net solution using the sparse ADMM algorithm
#' 
#' @param X sparse design matrix
#' @param y response vector
#' @param w L2 weight vector 
#' @param nu L1 weight vector 
#' @param D sparse penalty matrix
#' @param beta0 optional initial value for beta
#' @param z0 optional initial value for the primal variables
#' @param u0 optional initial value for the dual variables
#' @param lambda penalty parameter
#' @param alpha elastic-net mixing parameter 
#' @param eps (NOT USED) intercept specific penalty parameter
#' @param rho augmented Lagrangian penalty parameter
#' @param tau dual over-relaxation parameter
#' @param gamma proximal regularization parameter
#' @param intercept (NOT USED) if `TRUE`, do not penalize per first column of `X`
#' @param spthr proportion of zero entries to transform `X` and `D` into sparse matrices
#' @param smw (NOT USED) if `TRUE`, allows for the SMW solver in the `p>n` regime
#' @param precondition if `TRUE`, allows for diagonal preconditioning
#' @param objtol objective tolerance stopping criterion
#' @param reltol relative tolerance stopping criterion
#' @param abstol absolute tolerance stopping criterion
#' @param maxiter maximum number of iterations
#' @param verbose if `TRUE`, print the current optimization status
#' @param freq frequency of printing of the optimization status
#' 
#' @return 
#' \itemize{
#'   \item \code{beta}: estimated primal coefficients
#'   \item \code{z}: estimated primal variable
#'   \item \code{u}: estimated dual variables
#'   \item \code{r}: primal residuals
#'   \item \code{s}: dual residuals
#'   \item \code{niter}: number of iterations to reach convergence
#'   \item \code{trace}: data-frame collecting the optimization history
#' }
#' 
#' @keywords internal
solve_admm_spenet_PQ <- function(X, y, w, nu, D=diag(ncol(X)), d=numeric(nrow(D)),
                                 beta0=NULL, z0=NULL, u0=NULL, lambda=1.0, alpha=0.99, 
                                 eps=1e-8, rho=1.0, tau=0.0, gamma=0.0, intercept=TRUE, 
                                 spthr=0.9, smw=FALSE, precondition=FALSE,
                                 objtol=1e-5, reltol=1e-2, abstol=1e-4, 
                                 maxiter=1000, verbose=FALSE, freq=10) {
  
  # Model dimensions
  m <- nrow(D)
  n <- nrow(X)
  p <- ncol(X)
  
  sqrtm <- sqrt(m)
  sqrtn <- sqrt(n)
  sqrtp <- sqrt(p)
  sqrtnm <- sqrt(n+m)
  
  # Check for sparsity
  D = Matrix::Matrix(D, sparse=TRUE)
  X = Matrix::Matrix(X, sparse=TRUE)
  
  # Set the penalty vectors
  if ((alpha>1) || (alpha<=0)) {
    stop("ADMM-PQ-ENET: alpha must be real value in the interval (0,1].")
  }
  
  pL1 <- rep(lambda*alpha, m)
  pL2 <- rep(lambda*(1-alpha), m)
  
  # Preconditioning
  pwts <- c(rep(1/lambda, times=n), rep(1, times=m))
  if (precondition) {
    pwts[ (1:n)] <- 1 / (1e-6 + Matrix::rowSums(abs(X))*abs(nu)*lambda)
    pwts[-(1:n)] <- 1 / (1e-6 + Matrix::rowSums(abs(D)))
  }
  
  # Initialization
  # ... Primal variables: regression coefficients
  beta <- if (!is.null(beta0)) beta0 else rnorm(p, sd=0.1)
  # ... Primal variables: linear predictor
  eta <- as.vector(X %*% beta)
  NXb <- as.vector(nu * eta)
  Db <- as.vector(D %*% beta)
  Hb <- c(NXb, Db - d)
  # ... Primal variables: augmented variables
  z <- if (!is.null(z0)) z0 else Hb
  # ... Dual variables: Lagrange multipliers
  u <- if (!is.null(u0)) u0 else numeric(n+m) 
  # ... Primal variables: previous iteration
  z_old <- z
  # ... Primal and dual residuals
  r <- numeric(n+m)
  s <- numeric(n+m)
  
  # Pre-compute static quantities
  # ... Working weights, weighted responses, multiplicative factors
  wts <- w + rho*nu^2*pwts[1:n]
  d0 <- c(rep(0, n), d)
  wy_0 <- c(w * y, rep(0, m))
  rho_nu_0 <- rho * pwts * c(nu, rep(1, m))
  pen_0 <- c(rep(1/rho, n), pL1/rho) / pwts
  # ... System matrix factorization
  P <- Matrix::crossprod(D, (rho*pwts[-(1:n)]+pL2)*D) + diag(gamma,p,p)
  Q <- Matrix::crossprod(X, wts*X) + P
  cholQ <- Matrix::Cholesky(Q, pivot=TRUE)
  
  # Objective function 
  loss_PQ <- 0.5 * sum(w*(y - eta)^2) + sum(abs(NXb))
  pen_ridge <- 0.5 * sum(pL2*(Db - d)^2)
  pen_lasso <- sum(pL1*abs(Db - d))
  objval <- loss_PQ + pen_ridge + pen_lasso
  objval_old <- objval
  objval_min <- objval
  
  # Set the optimization history
  trace <- as.data.frame(matrix(NaN, nrow=maxiter, ncol=10))
  colnames(trace) <- c("iter", "objval", "diffobj", "norm_b", "norm_z", "norm_u", 
                       "norm_r", "norm_s", "eps_pri", "eps_dual")
  
  # Store the initial values
  trace$iter[1] <- 1
  trace$objval[1] <- objval
  trace$diffobj[1] <- NaN
  trace$norm_b[1] <- sqrt(sum(beta^2)) / sqrtp
  trace$norm_z[1] <- sqrt(sum(z^2)) / sqrtnm
  trace$norm_u[1] <- sqrt(sum(u^2)) / sqrtnm
  trace$norm_r[1] <- NaN
  trace$norm_s[1] <- NaN
  trace$eps_pri[1] <- NaN
  trace$eps_dual[1] <- NaN
  
  # Initial output
  iter = 1
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "\n")
  }
  
  # ADMM optimization loop
  for (iter in 2:maxiter) {
    
    # Primal update: beta
    e <- as.vector(Matrix::crossprod(rbind(X, D), wy_0 + rho_nu_0 * (d0+z-u))) + gamma*beta
    beta[] <- sp_chol_solve(cholQ, e)

    # Primal update: eta, NXb, Db
    eta[] <- as.vector(X %*% beta)
    NXb[] <- as.vector(nu * eta)
    Db[] <- as.vector(D %*% beta)
    Hb[] <- c(NXb, Db - d)
    
    # I dual update: u
    u[] <- u + tau * (Hb - z)
    
    # Primal update: z
    z_old[] <- z
    z[] <- soft_threshold(c(Hb + u), pen_0)
    
    # II dual update: u
    u[] <- u + (1-tau) * (Hb - z)
    
    # Primal and dual residuals
    r[] <- Hb - z
    s[] <- - rho * pwts * (z - z_old)
    
    # Convergence diagnostics
    norm_b <- sqrt(sum(beta^2))
    norm_z <- sqrt(sum(z^2))
    norm_u <- sqrt(sum(u^2))
    norm_r <- sqrt(sum(r^2))
    norm_s <- sqrt(sum(s^2))
    eps_pri <- max(sqrtp, sqrtnm) * abstol + reltol * max(norm_b, norm_z)
    eps_dual <- sqrtnm * abstol + reltol* rho * norm_u
    
    # Objective function
    loss_PQ <- 0.5 * sum(w*(y - eta)^2) + sum(abs(NXb))
    pen_ridge <- 0.5 * sum(pL2*(Db - d)^2)
    pen_lasso <- sum(pL1*abs(Db - d))
    objval <- loss_PQ + pen_ridge + pen_lasso
    diffobj <- abs(objval-objval_old) / (abs(objval)+1e-8)
    objval_min <- min(objval_min, objval_old)
    objval_old <- objval
    
    # Optimization state allocation
    trace$iter[iter] <- iter
    trace$objval[iter] <- objval
    trace$diffobj[iter] <- diffobj
    trace$norm_b[iter] <- norm_b / sqrtp
    trace$norm_z[iter] <- norm_z / sqrtm
    trace$norm_u[iter] <- norm_u / sqrtm
    trace$norm_r[iter] <- norm_r / sqrtm
    trace$norm_s[iter] <- norm_s / sqrtm
    trace$eps_pri[iter] <- eps_pri
    trace$eps_dual[iter] <- eps_dual
    
    # Intermediate output
    if (verbose && iter%%freq==0) {
      cat("  iter:", gettextf("%3d", iter),
          "  objval:", gettextf("%.4f", objval/n),
          "  diffobj:", gettextf("%.5f", diffobj),
          "  norm_r:", gettextf("%.4f", norm_r),
          "  norm_s:", gettextf("%.4f", norm_s),
          "\n")
    }
    
    # Convergence check
    chech_obj_sign <- (objval < objval_min)
    check_obj_diff <- (diffobj < objtol)
    check_r_norm <- (norm_r < eps_pri)
    check_s_norm <- (norm_s < eps_dual)
    if (chech_obj_sign && check_obj_diff && check_r_norm && check_s_norm) {
      break
    }
  }
  
  # Remove unused iterations
  if (iter < maxiter){
    trace <- trace[1:iter,]
  }
  
  # Return the results
  return(list(beta=beta, z=z, u=u, r=r, s=s, niter=iter, trace=trace))
}

#' @rdname solve_admm_spenet_PQ
#' @keywords internal
solve_admm_spenet_WLS <- function(X, y, w, D=diag(ncol(X)), d=numeric(nrow(D)),
                                  beta0=NULL, z0=NULL, u0=NULL, lambda=1.0, alpha=0.99, 
                                  eps=1e-8, rho=1.0, tau=0.0, gamma=0.0, intercept=TRUE, 
                                  spthr=0.9, smw=FALSE, precondition=FALSE,
                                  objtol=1e-5, reltol=1e-2, abstol=1e-4, 
                                  maxiter=1000, verbose=FALSE, freq=10) {
  
  # Model dimensions
  m <- nrow(D)
  n <- nrow(X)
  p <- ncol(X)
  
  sqrtm <- sqrt(m)
  sqrtn <- sqrt(n)
  sqrtp <- sqrt(p)
  
  # Check for sparsity
  D = Matrix::Matrix(D, sparse=TRUE)
  X = Matrix::Matrix(X, sparse=TRUE)
  
  # Set the penalty vectors
  if ((alpha>1) || (alpha<=0)) {
    stop("ADMM-WLS-ENET: alpha must be real value in the interval (0,1].")
  }
  
  pL1 <- rep(lambda*alpha, m)
  pL2 <- rep(lambda*(1-alpha), m)
  
  # Preconditioning
  if (precondition) {
    pwts <- 1 / (Matrix::rowSums(abs(D)) + 1e-6)
  } else {
    pwts <- rep(1, times=m)
  }
  
  # Initialization
  # ... Primal variables: regression coefficients
  beta <- if (!is.null(beta0)) beta0 else rnorm(p, sd=0.1)
  # ... Primal variables: linear predictor
  eta <- as.vector(X %*% beta)
  Db <- as.vector(D %*% beta)
  # ... Primal variables: augmented variables
  z <- if (!is.null(z0)) z0 else (Db - d)
  # ... Dual variables: Lagrange multipliers
  u <- if (!is.null(u0)) u0 else numeric(m)
  # ... Primal variables: previous iteration
  z_old <- z
  # ... Primal and dual residuals
  r <- numeric(m)
  s <- numeric(m)
  
  # Pre-compute static quantities
  # ... Working weights, weighted responses, multiplicative factors
  XWy <- as.vector(Matrix::crossprod(X, w * y))
  # ... System matrix factorization
  P <- Matrix::crossprod(sqrt(rho*pwts+pL2)*D) + diag(gamma,p,p)
  Q <- Matrix::crossprod(sqrt(w)*X) + P
  cholQ <- Matrix::Cholesky(Q, perm=TRUE)
  
  
  # Objective function 
  loss_WLS <- 0.5 * sum(w*(y - eta)^2)
  pen_ridge <- 0.5 * sum(pL2*(Db - d)^2)
  pen_lasso <- sum(pL1*abs(Db - d))
  objval <- loss_WLS + pen_ridge + pen_lasso
  objval_old <- objval
  objval_min <- objval
  
  # Set the optimization history
  trace <- as.data.frame(matrix(NaN, nrow=maxiter, ncol=10))
  colnames(trace) <- c("iter", "objval", "diffobj", "norm_b", 
                       "norm_z", "norm_u", "norm_r", "norm_s", 
                       "eps_pri", "eps_dual")
  
  # Store the initial values
  trace$iter[1] <- 1
  trace$objval[1] <- objval
  trace$diffobj[1] <- NaN
  trace$norm_b[1] <- sqrt(sum(beta^2)) / sqrtp
  trace$norm_z[1] <- sqrt(sum(z^2)) / sqrtm
  trace$norm_u[1] <- sqrt(sum(u^2)) / sqrtm
  trace$norm_r[1] <- NaN
  trace$norm_s[1] <- NaN
  trace$eps_pri[1] <- NaN
  trace$eps_dual[1] <- NaN
  
  # Initial output
  iter = 1
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "\n")
  }
  
  # ADMM optimization loop
  for (iter in 2:maxiter) {
    
    # Primal update: beta
    e <- XWy + rho * as.vector(Matrix::crossprod(D, pwts*(d+z-u))) + gamma*beta
    beta[] <- sp_chol_solve(cholQ, e)
    
    # Primal update: eta, NXb, Db
    eta[] <- as.vector(X %*% beta)
    Db[] <- as.vector(D %*% beta)
    
    # I dual update: u 
    u[] <- u + tau * (Db - d - z)
    
    # Primal update: z
    z_old[] <- z
    z[] <- soft_threshold(c(Db - d + u), (lambda*alpha/rho)/pwts)
    
    # II dual update: u
    u[] <- u + (1-tau) * (Db - d - z)
    
    # Primal and dual residuals
    r[] <- Db - d - z
    s[] <- - rho * pwts * (z - z_old)
    
    # Convergence diagnostics
    norm_b <- sqrt(sum(beta^2))
    norm_z <- sqrt(sum(z^2))
    norm_u <- sqrt(sum(u^2))
    norm_r <- sqrt(sum(r^2))
    norm_s <- sqrt(sum(s^2))
    eps_pri <- abstol * max(sqrtp, sqrtm) + reltol * max(norm_b, norm_z)
    eps_dual <- abstol * sqrtm + reltol * rho * norm_u
    
    # Objective function
    loss_WLS <- 0.5 * sum(w*(y - eta)^2)
    pen_ridge <- 0.5 * sum(pL2*(Db - d)^2)
    pen_lasso <- sum(pL1*abs(Db - d))
    objval <- loss_WLS + pen_ridge + pen_lasso
    diffobj <- abs(objval-objval_old) / (abs(objval)+1e-8)
    objval_min <- min(objval_min, objval_old)
    objval_old <- objval
    
    # Optimization state allocation
    trace$iter[iter] <- iter
    trace$objval[iter] <- objval
    trace$diffobj[iter] <- diffobj
    trace$norm_b[iter] <- norm_b / sqrtp
    trace$norm_z[iter] <- norm_z / sqrtm
    trace$norm_u[iter] <- norm_u / sqrtm
    trace$norm_r[iter] <- norm_r / sqrtm
    trace$norm_s[iter] <- norm_s / sqrtm
    trace$eps_pri[iter] <- eps_pri
    trace$eps_dual[iter] <- eps_dual
    
    # Intermediate output
    if (verbose && iter%%freq==0) {
      cat("  iter:", gettextf("%3d", iter),
          "  objval:", gettextf("%.4f", objval/n),
          "  diffobj:", gettextf("%.5f", diffobj),
          "  norm_r:", gettextf("%.4f", norm_r),
          "  norm_s:", gettextf("%.4f", norm_s),
          "\n")
    }
    
    # Convergence check
    chech_obj_sign <- (objval < objval_min)
    check_obj_diff <- (diffobj < objtol)
    check_r_norm <- (norm_r < eps_pri)
    check_s_norm <- (norm_s < eps_dual)
    if (chech_obj_sign && check_obj_diff && check_r_norm && check_s_norm) {
      break
    }
  }
  
  # Intermediate output
  if (verbose) {
    cat("  iter:", gettextf("%3d", iter),
        "  objval:", gettextf("%.4f", objval/n),
        "  diffobj:", gettextf("%.5f", diffobj),
        "  norm_r:", gettextf("%.4f", norm_r),
        "  norm_s:", gettextf("%.4f", norm_s),
        "\n")
  }
  
  # Remove unused iterations
  if (iter < maxiter){
    trace <- trace[1:iter,]
  }
  
  # Return the results
  return(list(beta=beta, z=z, u=u, r=r, s=s, niter=iter, trace=trace))
}

