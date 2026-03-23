
#' @title Ridge logistic regression via PG Gibbs sampling
#'
#' @param y A \verb{n} dimensional binary vector
#' @param X A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param D A \verb{m x p} dimensional penalty matrix.
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real) Regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param solver (string) Solver used for the linear systems and implicit matrix inversion. Must be one of 'auto', 'chol', 'smw', or 'sparse'
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param burn (integer) Burn-in period.
#' @param thin (integer) Thinning frequency.
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_mcmc_adj <- function(y, X, D=diag(ncol(X)), beta_start=NULL,
                               lambda=1.0, eps=1e-10, intercept=FALSE,
                               hpar=list(A=1.0, B=1.0, v=3),
                               prior=c("auto", "IG", "HW"), 
                               solver=c("auto", "chol", "smw", "sparse"),
                               scaling=TRUE, nscaling=10,
                               maxiter=5000L, burn=floor(maxiter/2), thin=5L,
                               verbose=FALSE, freq=floor(maxiter/20)){
  
  # Package for PG-RND
  library(BayesLogit)
  
  # Check the prior and solver
  prior <- match.arg(prior)
  solver <- match.arg(solver)
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the RND solver
  if (solver=="auto") {
    check_sp <- is(X, "sparseMatrix")
    check_np <- (p > n)
    
    if (check_sp) solver <- "sparse"
    if (!check_sp & !check_np) solver <- "chol"
    if (!check_sp &  check_np) solver <- "smw"
  }
  
  # Set the penalty parameters
  if (any(lambda<=0))
    stop("MCMC-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  lambda0 <- lambda
  pL2 <- rep(lambda, m)
  if (intercept) pL2[1] <- eps

  # Memory allocation
  trace <- list(
    beta = matrix(NA, nrow=maxiter, ncol=p),
    lambda = rep(NA, times=maxiter),
    alpha = rep(NA, times=maxiter),
    loglik = rep(NA, times=maxiter),
    logprior = rep(NA, times=maxiter),
    logpost = rep(NA, times=maxiter))
  
  summary <- list(
    lambda = list(mean=NA, inv=NA),
    beta = list(mean=rep(NA, p), var=rep(NA, p)),
    eta = list(mean=rep(NA, n), var=rep(NA, n)))
  
  setting = list(
    lambda = lambda,
    eps = eps, 
    intercept = intercept,
    hpar = hpar,
    prior = prior,
    solver = solver,
    scaling = scaling,
    nscaling = nscaling,
    maxiter = maxiter,
    burn = burn,
    thin = thin
  )
  
  # Parameter initialization
  iter <- 1
  
  sigma <- 0.3
  gamma <- 0.1
  avgrate <- 0.44
  
  # ... Regression parameters
  if (!is.null(beta_start)) {
    beta <- beta_start
  } else {
    if (solver=="chol") {
      Xty <- as.vector(crossprod(X, y-0.5))
      R <- chol(crossprod(rbind(0.5*X, sqrt(pL2)*D)))
      beta <- chol_rnd_sim(R, Xty)
    }
    if (solver=="smw") {
      beta <- smw_rnd_sim(0.5*X, (y-0.5)/0.5, pL2)
    }
    if (solver=="sparse") {
      Xty <- as.vector(Matrix::crossprod(X, y-0.5))
      Q <- Matrix::crossprod(rbind(0.5*X, sqrt(pL2)*D))
      R <- Matrix::Cholesky(Q, perm=FALSE, LDL=FALSE)
      beta <- sp_chol_rnd_sim(R, Xty)
    }
  }
  
  # ... Linear predictor
  eta <- as.vector(X %*% beta)
  delta <- as.vector(D %*% beta)
  
  # ... PG augmented variables
  w <- BayesLogit::rpg(n, h=1, z=eta)
  
  # ... Penalty parameter
  A_post <- hpar$A + .5*m
  B_post <- hpar$B + .5*sum(delta^2)
  lambda <- rgamma(1, shape=A_post, rate=B_post)
  
  if ( intercept) pL2[2:m] <- lambda
  if (!intercept) pL2[1:m] <- lambda
  
  # ... Log-posterior
  loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
  logpr <- 0.5 * (sum(log(pL2)) - m*log(2*pi) - sum(pL2*delta^2)) 
  logpr <- logpr + dgamma(lambda, shape=hpar$A, rate=hpar$B, log=TRUE)
  logpost <- loglik + logpr
  
  # Store the initial values
  trace$beta[iter,] <- beta
  trace$lambda[iter] <- lambda
  trace$loglik[iter] <- loglik
  trace$logprior[iter] <- logpr
  trace$logpost[iter] <- logpost
  
  # Progress bar
  if (verbose) {
    todo <- paste0(rep(".", floor(maxiter/freq)), collapse="")
    cat(paste0("   0/100 |", todo, "|", collapse=""))
  }
  
  # Gibbs sampling cycle
  for(iter in 2:maxiter) {
    # Progress bar
    if (verbose && iter%%freq==0) {
      perc <- gettextf("\r %3.0d/100", floor(100*iter/maxiter))
      done <- paste0(rep("=", floor(iter/freq)), collapse="")
      todo <- paste0(rep(".", floor((maxiter-iter)/freq)), collapse="")
      cat(paste0(perc, " |", done, todo, "|", collapse=""))
    }
    
    # Sample the regression parameters
    if (solver=="chol") {
      R <- chol(crossprod(rbind(sqrt(w)*X, sqrt(pL2)*D)))
      beta <- chol_rnd_sim(R, Xty)
    } 
    if (solver=="smw") {
      beta <- smw_rnd_sim(sqrt(w)*X, (y-0.5)/sqrt(w), pL2)
    }
    if (solver=="sparse") {
      Q <- Matrix::crossprod(rbind(sqrt(w)*X, sqrt(pL2)*D))
      R <- Matrix::Cholesky(Q, perm=FALSE, LDL=FALSE)
      beta <- sp_chol_rnd_sim(R, Xty)
    }
    
    # Update the linear predictor
    eta <- as.vector(X %*% beta)
    delta <- as.vector(D %*% beta)
    
    # Sample the PG augmented variables
    w <- BayesLogit::rpg(n, h=1, z=eta)
    
    # Sample the penalty parameter
    B_post <- hpar$B + .5*sum(delta^2)
    lambda <- rgamma(1, shape=A_post, rate=B_post)
    
    if ( intercept) pL2[2:m] <- lambda
    if (!intercept) pL2[1:m] <- lambda
    
    if (scaling) {
      
      # gamma <- .1 / (1 + iter/50)^.6
      # logsigma <- log(sigma) + gamma * (avgrate - 0.44)
      # logsigma <- ifelse(logsigma<0, logsigma, 0)
      # sigma <- exp(logsigma)
      
      for (k in 1:nscaling) {
        beta_null <- mean(beta)
        beta_ort <- beta - beta_null
        
        alpha <- rnorm(1, mean=0, sd=sigma)
        
        lambda_new <- lambda * exp(alpha)
        beta_new <- beta_null + beta_ort * exp(-alpha/2)
        
        eta_new <- as.vector(X %*% beta_new)
        delta_new <- as.vector(D %*% beta_new)
        
        loglik_old <- sum(y*eta - Rmpfr::log1pexp(eta))
        loglik_new <- sum(y*eta_new - Rmpfr::log1pexp(eta_new))
        logpr_rate <- hpar$A*alpha - hpar$B*lambda*(exp(alpha)-1)
        logacc_rate <- loglik_new - loglik_old + logpr_rate
        
        accrate <- min(1,exp(logacc_rate))
        if (k == 1) avgrate <- accrate
        if (k >  1) avgrate <- (1-1/k)*avgrate + (1/k)*accrate
        
        if (runif(1) < accrate) {
          lambda <- lambda_new
          beta <- beta_new
          eta <- eta_new
          delta <- delta_new
        }
      }
    }
    
    # Update log-posterior
    loglik <- sum(y*eta - Rmpfr::log1pexp(eta))
    logpr <- 0.5 * (sum(log(pL2)) - m*log(2*pi) - sum(pL2*delta^2)) 
    logpr <- logpr + dgamma(lambda, shape=hpar$A, rate=hpar$B, log=TRUE)
    logpost <- loglik + logpr
    
    # Store the result
    trace$beta[iter,] <- beta
    trace$lambda[iter] <- lambda
    trace$alpha[iter] <- avgrate
    trace$loglik[iter] <- loglik
    trace$logprior[iter] <- logpr
    trace$logpost[iter] <- logpost
  }
  
  if (verbose) cat("\n")
  
  # Posterior summaries
  idx <- seq(from=burn+1, to=maxiter, by=thin)
  
  summary$lambda$mean <- mean(trace$lambda[idx])
  summary$lambda$inv  <- mean(1/trace$lambda[idx])
  summary$beta$mean <- as.vector(apply(trace$beta[idx,], 2, mean))
  summary$beta$var  <- as.vector(apply(trace$beta[idx,], 2, var))
  summary$eta$mean <- as.vector(apply(X %*% t(trace$beta[idx,]), 1, mean))
  summary$eta$var  <- as.vector(apply(X %*% t(trace$beta[idx,]), 1, var))
  
  # Output
  return(list(trace=trace, summary=summary, call=setting))
}




#' @title Ridge logistic regression via variational approximation
#'
#' @param y A \verb{n} dimensional binary vector
#' @param X A \verb{n x p} dimensional design matrix. The intercept must be manually included
#' @param D A \verb{m x p} dimensional penalty matrix.
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients
#' @param lambda (real) Regularization parameters for ridge regression
#' @param eps (real) Small value added to the intercept term to avoid numerical issues
#' @param intercept (boolean) Is the intercept included in `X` or not?
#' @param solver (string) Solver used for the linear systems and implicit matrix inversion. Must be one of 'auto', 'chol', 'smw', or 'sparse'
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error
#' @param abstol (real) Convergence threshold for the absolute change in the log-likelihood
#' @param reltol (real) Convergence threshold for the relative change in the log-likelihood
#' @param verbose (boolean) Print the intermediate state of the optimization
#' @param freq (int) How often print the optimization state
#' 
#' @return The function returns a list with optimal coefficients and likelihood path through MM iterations
#' 
#' @export
fit_logit_mfvb_adj <- function(y, X, D=diag(ncol(X)), type=c("PQ", "PG", "BL"), 
                               beta_start=NULL, lambda=1.0, eps=1e-10, intercept=FALSE,  
                               hpar=list(A=1.0, B=1.0, v=3),
                               prior=c("auto", "IG", "HW"), 
                               solver=c("auto", "chol", "smw", "sparse"),
                               maxiter=1000L, abstol=1e-4, reltol=1e-4,
                               verbose=FALSE, freq=10L){
  
  # Set the VB method
  type <- match.arg(type)
  prior <- match.arg(prior)
  
  if (!is.null(D) && solver=="smw")
    stop("Method 'smw' does not allow for a general penalty matrix D.")
  
  # Set the fitting function
  fit <- switch(type,
                "PQ" = fit_logit_pqvb_adj, 
                "PG" = fit_logit_pgvb_adj, 
                "BL" = fit_logit_blvb_adj)
  
  # Run the optimization
  fit(y=y, X=X, D=D, beta_start=beta_start, lambda=lambda, eps=eps, 
      intercept=intercept, hpar=hpar, prior=prior, solver=solver, maxiter=maxiter, 
      abstol=abstol, reltol=reltol, verbose=verbose, freq=freq)
}


#' @rdname fit_logit_mfvb
#' @export
fit_logit_pqvb_adj <- function(y, X, D=diag(ncol(X)), beta_start=NULL, 
                               lambda=1.0, eps=1e-10, intercept=FALSE, 
                               hpar=list(A=1.0, B=1.0, v=3),
                               prior=c("auto", "IG", "HW"), 
                               solver=c("auto", "chol", "smw", "sparse"),
                               maxiter=1000L, abstol=1e-4, reltol=1e-4,
                               verbose=FALSE, freq=10L){
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the RND solver
  solver <- match.arg(solver)
  if (solver=="auto") {
    check_sp <- is(X, "sparseMatrix")
    check_np <- (p > n)
    
    if (check_sp) solver <- "sparse"
    if (!check_sp & !check_np) solver <- "chol"
    if (!check_sp &  check_np) solver <- "smw"
  }
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("PQVB-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda, m)
  if (intercept) pL2[1] <- eps
  
  # Parameter initialization
  mu <- beta_start
  if (solver=="chol") {
    # ... Regression parameters
    R <- chol(crossprod(rbind(0.5*X, sqrt(lambda)*D)))
    invR <- solve(R)
    if (is.null(beta_start))
      mu <- chol_solve(R, crossprod(X, 2*(y-0.5)))
    
    # ... Squared q-expectation of beta
    trVP <- sum((D %*% invR)^2)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -2*sum(log(diag(R)))
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    sq_eta <- sqrt(rowSums((X %*% invR)^2))
  } 
  if (solver=="smw") {
    # ... VB weights and pseudo-data
    w_VB <- rep(0.25, times=n)
    z_VB <- 2*(y-0.5)
    
    # ... Sufficient statistics
    XPXt <- X %*% (t(X) / pL2)
    XtWz <- as.vector(crossprod(X, z_VB))
    
    # ... Regression parameters
    R <- chol(diag(1/w_VB) + XPXt)
    if (is.null(beta_start))
      mu <- smw_chol_solve(R, X, pL2, XtWz)
    
    # ... Squared q-expectation of beta
    Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
    Vq_logdet <- smw_chol_logdet(R, w_VB, pL2)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    sq_eta <- sqrt(colSums(t(X) * smw_chol_solve(R, X, pL2, t(X))))
  }
  if (solver=="sparse") {
    # ... Regression parameters
    P <- Matrix::crossprod(D)
    Q <- Matrix::crossprod(0.5*X) + lambda*P
    R <- Matrix::Cholesky(Q, perm=TRUE)
    XtWz <- as.vector(Matrix::crossprod(X, 2*(y-0.5)))
    if (is.null(beta_start))
      mu <- sp_chol_solve(R, XtWz)
    
    # ... Squared q-expectation of beta
    trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -sp_chol_logdet(R)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    sq_eta <- sqrt(Matrix::colSums(t(X) * sp_chol_solve(R, t(X))))
  }
  
  # ... Penalty parameter
  Aq_lambda <- hpar$A + .5*m
  Bq_lambda <- hpar$B + .5*Eq_sq_beta
  Eq_lambda <- Aq_lambda / Bq_lambda
  Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
  Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
  
  if ( intercept) pL2[2:m] <- Eq_lambda
  if (!intercept) pL2[1:m] <- Eq_lambda
  
  # ... Gaussian pdf and cdf of eta in 0
  pdf_eta <- dnorm(mq_eta, mean=0, sd=sq_eta)
  cdf_eta <- pnorm(mq_eta, mean=0, sd=sq_eta)
  
  # ... Squared and absolute q-expectation of eta
  Eq_sq_eta <- mq_eta*mq_eta + sq_eta*sq_eta
  Eq_abs_eta <- mq_eta*(2*cdf_eta-1) + 2*(sq_eta*sq_eta)*pdf_eta
  Eq_sgn_eta <- 2*cdf_eta - 1
  Eq_dlt_eta <- 2*pdf_eta
  
  # ... PQ local variables
  xi <- Eq_sq_eta / Eq_abs_eta
  
  # ... L2 and L1 PQ weights
  w_PG <- tanh(xi/2) / (2*xi)
  w_PQ <- 2*w_PG - 2*log(cosh(xi/2)) / (xi*xi)
  v_PQ <- abs(xi) * (w_PG - w_PQ)
  h_PQ <- 0.5*xi - Rmpfr::log1pexp(xi)
  h_VB <- sum(h_PQ - 0.5*w_PQ*(Eq_sq_eta-xi*xi) - v_PQ*(Eq_abs_eta-xi))
  
  # ... VB log-likelihood
  Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
  
  # ... ELBO
  # Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
  # Eq_logp <- -0.5*(Eq_sq_beta - p*sum(log(pL2)))
  # Eq_logq <- -0.5*Vq_logdet
  # objval <- Eq_loglik + Eq_logp - Eq_logq
  objval <- Eq_loglik + 0.5*Vq_logdet
  objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
  objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)
  
  # Set the initial optimization state
  iter <- 1
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_mu <- NaN
  diff_rel_eta <- NaN
  diff_rel_xi <- NaN
  
  # Set the optimization control parameters
  setting <- list(
    type = "PQ", lambda = lambda, eps = eps, intercept = intercept,
    hpar = hpar, prior = prior, solver = solver, 
    maxiter = maxiter, abstol = abstol, reltol = reltol)
  
  # Set the optimization summary
  summary <- list(beta=list(), eta=list(), lambda=list(), xi=NULL, elbo=NULL)
  state <- list(xi=NULL, w_PQ=NULL, v_PQ=NULL, h_PQ=NULL, h_VB=NULL, elbo=NULL, niter=NULL)
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=8))
  colnames(trace) <- c("iter", "objval", "lambda", "diff_abs_obj", "diff_rel_obj", 
                       "diff_rel_mu", "diff_rel_eta", "diff_rel_xi")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                    diff_rel_mu, diff_rel_eta, diff_rel_xi)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "mu"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "xi"), "\n",
        paste0(rep("-", 73), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        "\n", sep="")
  }
  
  # Iterative procedure
  for(iter in 2:maxiter){
    
    ## Set the old parameters
    objval_old <- objval
    eta_old <- mq_eta
    mu_old <- mu
    xi_old <- xi
    
    # VB weights and pseudo-data
    w_VB <- w_PQ + v_PQ*Eq_dlt_eta
    z_VB <- ((y-0.5) - w_PQ*mq_eta - v_PQ*Eq_sgn_eta)/w_VB + mq_eta
    
    # Update regression parameters: beta
    if (solver=="chol") {
      # ... Regression parameters
      R <- chol(crossprod(rbind(sqrt(w_VB)*X, sqrt(Eq_lambda)*D)))
      invR <- solve(R)
      XtWz <- crossprod(X, w_VB*z_VB)
      mu <- chol_solve(R, XtWz)
      
      # ... Squared q-expectation of beta
      trVP <- sum((D %*% invR)^2)
      Eq_sq_beta <- sum((D %*% mu)^2) + trVP
      Vq_logdet <- -2*sum(log(diag(R)))
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      sq_eta <- sqrt(rowSums((X %*% invR)^2))
    }
    if (solver=="smw") {
      # ... Regression parameters
      R <- chol(diag(1/w_VB) + XPXt)
      XtWz <- as.vector(crossprod(X, w_VB*z_VB))
      mu <- smw_chol_solve(R, X, pL2, XtWz)
      
      # ... Squared q-expectation of beta
      Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
      Vq_logdet <- smw_chol_logdet(R, w_VB, pL2)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      sq_eta <- sqrt(colSums(t(X) * smw_chol_solve(R, X, pL2, t(X))))
    }
    if (solver=="sparse") {
      # ... Regression parameters
      Q <- Matrix::crossprod(X, w_VB*X) + Eq_lambda*P
      R <- Matrix::Cholesky(Q, perm=TRUE)
      XtWz <- as.vector(Matrix::crossprod(X, w_VB*z_VB))
      mu <- sp_chol_solve(R, XtWz)
      
      # ... Squared q-expectation of beta
      trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
      Eq_sq_beta <- sum((D %*% mu)^2) + trVP
      Vq_logdet <- -sp_chol_logdet(R)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      sq_eta <- sqrt(Matrix::colSums(t(X) * sp_chol_solve(R, t(X))))
    }
    
    # ... Penalty parameter
    Aq_lambda <- hpar$A + .5*m
    Bq_lambda <- hpar$B + .5*Eq_sq_beta
    Eq_lambda <- Aq_lambda / Bq_lambda
    Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
    Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
    
    if ( intercept) pL2[2:m] <- Eq_lambda
    if (!intercept) pL2[1:m] <- Eq_lambda
    
    # ... Gaussian pdf and cdf of eta in 0
    pdf_eta <- dnorm(mq_eta, mean=0, sd=sq_eta)
    cdf_eta <- pnorm(mq_eta, mean=0, sd=sq_eta)
    
    # ... Squared and absolute q-expectation of eta
    Eq_sq_eta <- mq_eta*mq_eta + sq_eta*sq_eta
    Eq_abs_eta <- mq_eta*(2*cdf_eta-1) + 2*(sq_eta*sq_eta)*pdf_eta
    Eq_sgn_eta <- 2*cdf_eta - 1
    Eq_dlt_eta <- 2*pdf_eta
    
    # Update PQ local variables: xi
    xi <- Eq_sq_eta / Eq_abs_eta
    
    # Update L2 and L1 PQ weights
    w_PG <- tanh(xi/2) / (2*xi)
    w_PQ <- 2*w_PG - 2*log(cosh(xi/2)) / (xi*xi)
    v_PQ <- abs(xi) * (w_PG - w_PQ)
    h_PQ <- 0.5*xi - Rmpfr::log1pexp(xi)
    h_VB <- sum(h_PQ - 0.5*w_PQ*(Eq_sq_eta-xi*xi) - v_PQ*(Eq_abs_eta-xi))
    
    # ... VB log-likelihood
    Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
    
    # ... ELBO
    # Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
    # Eq_logp <- -0.5*(Eq_sq_beta - p*sum(log(pL2)))
    # Eq_logq <- -0.5*Vq_logdet
    # objval <- Eq_loglik + Eq_logp - Eq_logq
    objval <- Eq_loglik + 0.5*Vq_logdet
    objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
    objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)
    
    # Check convergence
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_mu  <- maxreldiff(mu, mu_old)
    diff_rel_eta <- maxreldiff(mq_eta, eta_old)
    diff_rel_xi  <- maxreldiff(xi, xi_old)
    
    # Store the current state
    trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                      diff_rel_mu, diff_rel_eta, diff_rel_xi)
    
    # Print the current optimization state
    if(verbose && iter%%freq==0){
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/(n+p)),
          gettextf(" %10.3e", diff_abs_obj/(n+p)),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_mu),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf(" %10.3e", diff_rel_xi),
          "\n", sep="")
    }
    
    if((diff_abs_obj<(n+p)*abstol) && (diff_rel_obj<reltol)){
      break
    }
  }
  
  # Compute the final omega
  omega = switch(solver,
                 "chol" = tcrossprod(invR),
                 "smw" = smw_chol_solve(R, X, pL2, diag(p)),
                 "sparse" = list(var=sparseinv::Takahashi_Davis(Q), cholQ=R))
  
  # Print the final optimization state
  if(verbose){
    cat(gettextf(" %5d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        gettextf(" %10.3e", diff_abs_obj/(n+p)),
        gettextf(" %10.3e", diff_rel_obj),
        gettextf(" %10.3e", diff_rel_mu),
        gettextf(" %10.3e", diff_rel_eta),
        gettextf(" %10.3e", diff_rel_xi), "\n",
        paste0(rep("-", 73), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("PQVB-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  # Set the final output
  summary$beta$mean <- mu
  summary$beta$var <- omega
  summary$beta$chol <- R
  summary$beta$sqtr <- Eq_sq_beta
  summary$beta$logdet <- Vq_logdet
  
  summary$eta$mean <- mq_eta
  summary$eta$var <- sq_eta^2
  
  summary$lambda$shape <- Aq_lambda
  summary$lambda$rate <- Bq_lambda
  summary$lambda$mean <- Eq_lambda
  summary$lambda$inv <- Eq_inv_lambda
  summary$lambda$log <- Eq_log_lambda
  
  state$xi <- xi
  state$w_PQ <- w_PQ
  state$v_PQ <- v_PQ
  state$h_PQ <- h_PQ
  state$h_VB <- h_VB
  state$elbo <- objval
  state$niter <- iter
  
  return(list(summary=summary, state=state, trace=trace[1:iter,], call=setting))
}

#' @rdname fit_logit_mfvb
#' @export
fit_logit_pgvb_adj <- function(y, X, D=diag(ncol(X)), beta_start=NULL, 
                               lambda=1.0, eps=1e-10, intercept=FALSE,
                               hpar=list(A=1.0, B=1.0, v=3),
                               prior=c("auto", "IG", "HW"), 
                               solver=c("auto", "chol", "smw", "sparse"),
                               maxiter=1000L, abstol=1e-4, reltol=1e-4,
                               verbose=FALSE, freq=10L){
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the RND solver
  solver <- match.arg(solver)
  if (solver=="auto") {
    check_sp <- is(X, "sparseMatrix")
    check_np <- (p > n)
    
    if (check_sp) solver <- "sparse"
    if (!check_sp & !check_np) solver <- "chol"
    if (!check_sp &  check_np) solver <- "smw"
  }
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("PGVB-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda, m)
  if (intercept) pL2[1] <- eps
  
  # Parameter initialization
  mu <- beta_start
  if (solver=="chol") {
    # ... Sufficient statistics
    Xty <- as.vector(crossprod(X, y-0.5))
    
    # ... Regression parameters
    R <- chol(crossprod(rbind(0.5*X, lambda*D)))
    invR <- solve(R)
    if (is.null(beta_start))
      mu <- chol_solve(R, 2*Xty)
    
    # ... Squared q-expectation of beta
    trVP <- sum((D %*% invR)^2)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -2*sum(log(diag(R)))
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- rowSums((X %*% invR)^2)
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  if (solver=="smw") {
    # ... Sufficient statistics
    XPXt <- X %*% (t(X) / pL2)
    Xty <- as.vector(crossprod(X, y-0.5))
    
    # ... Regression parameters
    R <- chol(diag(4.0, n, n) + XPXt)
    if (is.null(beta_start))
      mu <- smw_chol_solve(R, X, pL2, 2*Xty)
    
    # ... Squared q-expectation of beta
    Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
    Vq_logdet <- smw_chol_logdet(R, rep(0.25, n), pL2)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- colSums(t(X) * smw_chol_solve(R, X, pL2, t(X)))
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  if (solver=="sparse") {
    # ... Regression parameters
    # P <- Matrix::crossprod(sqrt(pL2)*D)
    P <- Matrix::crossprod(D)
    Q <- Matrix::crossprod(0.5*X) + lambda*P
    R <- Matrix::Cholesky(Q, perm=TRUE)
    Xty <- as.vector(Matrix::crossprod(X, y-0.5))
    if (is.null(beta_start))
      mu <- sp_chol_solve(R, 2*Xty)
    
    # ... Squared q-expectation of beta
    trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -sp_chol_logdet(R)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- Matrix::colSums(t(X) * sp_chol_solve(R, t(X)))
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  
  # ... Penalty parameter
  Aq_lambda <- hpar$A + .5*m
  Bq_lambda <- hpar$B + .5*Eq_sq_beta
  Eq_lambda <- Aq_lambda / Bq_lambda
  Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
  Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
  
  if ( intercept) pL2[2:m] <- Eq_lambda
  if (!intercept) pL2[1:m] <- Eq_lambda
  
  # ... PG local variable
  xi <- sqrt(Eq_sq_eta)
  
  # ... PG weights
  w_PG <- tanh(xi/2) / (2*xi)
  h_PG <- 0.5*xi - Rmpfr::log1pexp(xi)
  h_VB <- sum(h_PG - 0.5*w_PG*(Eq_sq_eta-xi*xi))
  
  # h_PQ <- get_hfun_PQ(mq_eta, sqrt(Vq_eta))
  
  # ... VB log-likelihood
  Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
  
  # ... ELBO
  # Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
  # Eq_logp <- -0.5*(Eq_sq_beta - p*sum(log(pL2)))
  # Eq_logq <- -0.5*Vq_logdet
  objval <- Eq_loglik + 0.5*Vq_logdet
  objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
  objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)
    
  # Set the initial optimization state
  iter <- 1
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_mu <- NaN
  diff_rel_eta <- NaN
  diff_rel_xi <- NaN
  
  # Set the optimization control parameters
  setting <- list(
    type = "PG", lambda = lambda, eps = eps, intercept = intercept,
    hpar = hpar, prior = prior, solver = solver, 
    maxiter = maxiter, abstol = abstol, reltol = reltol)
  
  # Set the optimization summary
  summary <- list(beta=list(), eta=list(), lambda=list(), xi=NULL, elbo=NULL)
  state <- list(xi=NULL, w_PG=NULL, h_PG=NULL, h_VB=NULL, elbo=NULL, niter=NULL)
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=8))
  colnames(trace) <- c("iter", "objval", "lambda", "diff_abs_obj", "diff_rel_obj", 
                       "diff_rel_mu", "diff_rel_eta", "diff_rel_xi")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                    diff_rel_mu, diff_rel_eta, diff_rel_xi)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "mu"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "xi"), "\n",
        paste0(rep("-", 73), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        "\n", sep="")
  }
  
  # Iterative procedure
  for(iter in 2:maxiter){
    
    ## Set the old parameters
    objval_old <- objval
    eta_old <- mq_eta
    mu_old <- mu
    xi_old <- xi
    
    ## Update regression parameters: beta
    if (solver=="chol") {
      # ... Regression parameters
      R <- chol(crossprod(rbind(sqrt(w_PG)*X, sqrt(Eq_lambda)*D)))
      invR <- solve(R)
      mu <- chol_solve(R, Xty)
      
      # ... Squared q-expectation of beta
      Eq_sq_beta <- sum((D %*% mu)^2) + sum((D %*% invR)^2)
      Vq_logdet <- -2*sum(log(diag(R)))
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- rowSums((X %*% invR)^2)
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    if (solver=="smw") {
      # ... Regression parameters
      R <- chol(diag(1/w_PG) + XPXt)
      mu <- smw_chol_solve(R, X, pL2, Xty)
      
      # ... Squared q-expectation of beta
      Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
      Vq_logdet <- smw_chol_logdet(R, w_PG, pL2)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- colSums(t(X) * smw_chol_solve(R, X, pL2, t(X)))
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    if (solver=="sparse") {
      # ... Regression parameters
      Q <- Matrix::crossprod(X, w_PG*X) + Eq_lambda*P
      R <- Matrix::Cholesky(Q, perm=TRUE)
      mu <- sp_chol_solve(R, Xty)
      
      # ... Squared q-expectation of beta
      trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
      Eq_sq_beta <- sum((D %*% mu)^2) + trVP
      Vq_logdet <- -sp_chol_logdet(R)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- Matrix::colSums(t(X) * sp_chol_solve(R, t(X)))
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    
    # ... Penalty parameter
    Aq_lambda <- hpar$A + .5*m
    Bq_lambda <- hpar$B + .5*Eq_sq_beta
    Eq_lambda <- Aq_lambda / Bq_lambda
    Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
    Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
    
    if ( intercept) pL2[2:m] <- Eq_lambda
    if (!intercept) pL2[1:m] <- Eq_lambda
    
    # Update local variable: xi
    xi <- sqrt(Eq_sq_eta)
    
    # Update L2 and L1 PQ weights
    w_PG <- tanh(xi/2) / (2*xi)
    h_PG <- 0.5*xi - Rmpfr::log1pexp(xi)
    h_VB <- sum(h_PG - 0.5*w_PG*(Eq_sq_eta-xi*xi))
    
    # h_PQ <- get_hfun_PQ(mq_eta, sqrt(Vq_eta))
    
    # Update VB log-likelihood
    Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
    
    # Update ELBO
    # Eq_logp <- -0.5*(Eq_sq_beta - p*sum(log(pL2)))
    # Eq_logq <- -0.5*Vq_logdet
    # objval <- Eq_loglik + Eq_logp - Eq_logq
    objval <- Eq_loglik + 0.5*Vq_logdet
    objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
    objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)
    
    
    # Check convergence
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_mu  <- maxreldiff(mu, mu_old)
    diff_rel_eta <- maxreldiff(mq_eta, eta_old)
    diff_rel_xi  <- maxreldiff(xi, xi_old)
    
    # Store the current state
    trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                      diff_rel_mu, diff_rel_eta, diff_rel_xi)
    
    # Print the current optimization state
    if(verbose && iter%%freq==0){
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/(n+p)),
          gettextf(" %10.3e", diff_abs_obj/(n+p)),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_mu),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf(" %10.3e", diff_rel_xi),
          "\n", sep="")
    }
    
    if((diff_abs_obj<(n+p)*abstol) && (diff_rel_obj<reltol)){
      break
    }
  }
  
  # Compute the final omega
  omega = switch(solver,
                 "chol" = tcrossprod(invR),
                 "smw" = smw_chol_solve(R, X, pL2, diag(p)),
                 "sparse" = list(var=sparseinv::Takahashi_Davis(Q), cholQ=R))
  
  # Print the final optimization state
  if(verbose){
    cat(gettextf(" %5d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        gettextf(" %10.3e", diff_abs_obj/(n+p)),
        gettextf(" %10.3e", diff_rel_obj),
        gettextf(" %10.3e", diff_rel_mu),
        gettextf(" %10.3e", diff_rel_eta),
        gettextf(" %10.3e", diff_rel_xi), "\n",
        paste0(rep("-", 73), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("PGVB-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  # Set the final output
  summary$beta$mean <- mu
  summary$beta$var <- omega
  summary$beta$chol <- R
  summary$beta$sqtr <- Eq_sq_beta
  summary$beta$logdet <- Vq_logdet
  
  summary$eta$mean <- mq_eta
  summary$eta$var <- Vq_eta
  
  summary$lambda$shape <- Aq_lambda
  summary$lambda$rate <- Bq_lambda
  summary$lambda$mean <- Eq_lambda
  summary$lambda$inv <- Eq_inv_lambda
  summary$lambda$log <- Eq_log_lambda
  
  state$xi <- xi
  state$w_PG <- w_PG
  state$h_PG <- h_PG
  state$h_VB <- h_VB
  state$elbo <- objval
  state$niter <- iter
  
  return(list(summary=summary, state=state, trace=trace[1:iter,], call=setting))
}


#' @rdname fit_logit_mfvb
#' @export
fit_logit_blvb_adj <- function(y, X, D=diag(ncol(X)), beta_start=NULL, 
                               lambda=1.0, eps=1e-10, intercept=FALSE, 
                               hpar=list(A=1.0, B=1.0, v=3),
                               prior=c("auto", "IG", "HW"), 
                               solver=c("auto", "chol", "smw", "sparse"),
                               maxiter=1000L, abstol=1e-4, reltol=1e-4,
                               verbose=FALSE, freq=10L){
  
  # Set the model dimensions
  p <- ncol(X)
  n <- nrow(X)
  m <- nrow(D)
  
  # Set the RND solver
  solver <- match.arg(solver)
  if (solver=="auto") {
    check_sp <- is(X, "sparseMatrix")
    check_np <- (p > n)
    
    if (check_sp) solver <- "sparse"
    if (!check_sp & !check_np) solver <- "chol"
    if (!check_sp &  check_np) solver <- "smw"
  }
  
  # Set the penalty parameters
  if (lambda<=0)
    stop("BLVB-LOGIT: lambda must be a positive real value.", call.=FALSE)
  
  pL2 <- rep(lambda, m)
  if (intercept) pL2[1] <- eps
  
  # Parameter initialization
  mu <- beta_start
  if (solver=="chol") {
    # ... Sufficient statistics
    XtWz <- as.vector(crossprod(X, y-0.5))
    
    # ... Regression parameters
    R <- chol(crossprod(rbind(0.5*X, sqrt(lambda)*D)))
    invR <- solve(R)
    if (is.null(beta_start))
      mu <- chol_solve(R, 2*XtWz)
    
    # ... Squared q-expectation of beta
    trVP <- sum((D %*% invR)^2)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -2*sum(log(diag(R)))
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- rowSums((X %*% invR)^2)
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  if (solver=="smw") {
    # ... Sufficient statistics
    XPXt <- X %*% (t(X) / pL2)
    XtWz <- as.vector(crossprod(X, y-0.5))
    
    # ... Regression parameters
    R <- chol(diag(4.0, n, n) + XPXt)
    if (is.null(beta_start))
      mu <- smw_chol_solve(R, X, pL2, 2*Xty)
    
    # ... Squared q-expectation of beta
    Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
    Vq_logdet <- smw_chol_logdet(R, rep(0.25, n), pL2)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- colSums(t(X) * smw_chol_solve(R, X, pL2, t(X)))
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  if (solver=="sparse") {
    # ... Regression parameters
    P <- Matrix::crossprod(D)
    Q <- Matrix::crossprod(0.5*X) + lambda*P
    R <- Matrix::Cholesky(Q, perm=TRUE)
    XtWz <- as.vector(Matrix::crossprod(X, y-0.5))
    if (is.null(beta_start))
      mu <- sp_chol_solve(R, 2*XtWz)
    
    # ... Squared q-expectation of beta
    trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
    Eq_sq_beta <- sum((D %*% mu)^2) + trVP
    Vq_logdet <- -sp_chol_logdet(R)
    
    # ... Linear predictor
    mq_eta <- as.vector(X %*% mu)
    Vq_eta <- Matrix::colSums(t(X) * sp_chol_solve(R, t(X)))
    Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
  }
  
  # ... Penalty parameter
  Aq_lambda <- hpar$A + .5*m
  Bq_lambda <- hpar$B + .5*Eq_sq_beta
  Eq_lambda <- Aq_lambda / Bq_lambda
  Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
  Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
  
  if ( intercept) pL2[2:m] <- Eq_lambda
  if (!intercept) pL2[1:m] <- Eq_lambda
  
  # ... BL local variable
  xi <- mq_eta
  
  # ... BL weights
  w_BL <- rep(0.25, times=n)
  z_BL <- 0.25*xi + y - 1/(1+exp(-xi))
  h_BL <- 0.5*xi - Rmpfr::log1pexp(xi)
  g_BL <- 0.5 - 1/(1+exp(-xi))
  h_VB <- sum(h_BL + g_BL*(mq_eta-xi) - 0.5*w_BL*(Eq_sq_eta-2*xi*mq_eta+xi*xi))
  
  # ... VB log-likelihood
  Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
  
  # ... ELBO
  objval <- Eq_loglik + 0.5*Vq_logdet
  objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
  objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)
  
  # Set the initial optimization state
  iter <- 1
  objval_old <- objval
  diff_abs_obj <- NaN
  diff_rel_obj <- NaN
  diff_rel_mu <- NaN
  diff_rel_eta <- NaN
  diff_rel_xi <- NaN
  
  # Set the optimization control parameters
  setting <- list(
    type="BL", lambda = lambda, eps = eps, intercept = intercept,
    hpar = hpar, prior = prior, solver = solver, 
    maxiter = maxiter, abstol = abstol, reltol = reltol)
  
  # Set the optimization summary
  summary <- list(beta=list(), eta=list(), lambda=list(), xi=NULL, elbo=NULL)
  state <- list(xi=NULL, w_PG=NULL, h_PG=NULL, h_VB=NULL, elbo=NULL, niter=NULL)
  
  # Set the optimization trace
  trace <- as.data.frame(matrix(NA, nrow=maxiter, ncol=8))
  colnames(trace) <- c("iter", "objval", "lambda", "diff_abs_obj", "diff_rel_obj", 
                       "diff_rel_mu", "diff_rel_eta", "diff_rel_xi")
  
  # Store the current state
  trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                    diff_rel_mu, diff_rel_eta, diff_rel_xi)
  
  # Print the initial optimization state
  if (verbose) {
    cat(gettextf("%6s", "iter"), 
        gettextf("%12s", "objval"),
        gettextf("%11s", "absdiff"),
        gettextf("%11s", "reldiff"),
        gettextf("%11s", "mu"),
        gettextf("%11s", "eta"),
        gettextf("%10s", "xi"), "\n",
        paste0(rep("-", 73), collapse=""), "\n",
        gettextf("%6d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        "\n", sep="")
  }
  
  # Iterative procedure
  for(iter in 2:maxiter){
    
    ## Set the old parameters
    objval_old <- objval
    eta_old <- mq_eta
    mu_old <- mu
    xi_old <- xi
    
    ## Update regression parameters: beta
    if (solver=="chol") {
      # ... Regression parameters
      R <- chol(crossprod(rbind(sqrt(w_BL)*X, sqrt(Eq_lambda)*D)))
      invR <- solve(R)
      XtWz <- as.vector(crossprod(X, z_BL))
      mu <- chol_solve(R, XtWz)
      
      # ... Squared q-expectation of beta
      trVP <- sum((D %*% invR)^2)
      Eq_sq_beta <- sum((D %*% mu)^2) + trVP
      Vq_logdet <- -2*sum(log(diag(R)))
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- rowSums((X %*% invR)^2)
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    if (solver=="smw") {
      # ... Regression parameters
      XtWz <- as.vector(crossprod(X, z_BL))
      mu <- smw_chol_solve(R, X, pL2, XtWz)
      
      # ... Squared q-expectation of beta
      Eq_sq_beta <- sum(pL2*mu^2) + sum(diag(chol_solve(R, XPXt)))
      Vq_logdet <- smw_chol_logdet(R, w_BL, pL2)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- colSums(t(X) * smw_chol_solve(R, X, pL2, t(X)))
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    if (solver=="sparse") {
      # ... Regression parameters
      Q <- Matrix::crossprod(X, w_BL*X) + Eq_lambda*P
      R <- Matrix::Cholesky(Q, perm=TRUE)
      XtWz <- as.vector(Matrix::crossprod(X, z_BL))
      mu <- sp_chol_solve(R, XtWz)
      
      # ... Squared q-expectation of beta
      trVP <- sp_chol_edf(R, P, rank=10, nsample=100, seed=1234)
      Eq_sq_beta <- sum((D %*% mu)^2) + trVP
      Vq_logdet <- -sp_chol_logdet(R)
      
      # ... Linear predictor
      mq_eta <- as.vector(X %*% mu)
      Vq_eta <- Matrix::colSums(t(X) * sp_chol_solve(R, t(X)))
      Eq_sq_eta <- mq_eta*mq_eta + Vq_eta
    }
    
    # ... Penalty parameter
    Aq_lambda <- hpar$A + .5*m
    Bq_lambda <- hpar$B + .5*Eq_sq_beta
    Eq_lambda <- Aq_lambda / Bq_lambda
    Eq_log_lambda <- digamma(Aq_lambda) - log(Bq_lambda)
    Eq_inv_lambda <- Bq_lambda / (Aq_lambda-1)
    
    if ( intercept) pL2[2:m] <- Eq_lambda
    if (!intercept) pL2[1:m] <- Eq_lambda
    
    # Update local variable: xi
    xi <- mq_eta
    
    # Update L2 and L1 PQ weights
    z_BL <- 0.25*xi + y - 1/(1+exp(-xi))
    h_BL <- 0.5*xi - Rmpfr::log1pexp(xi)
    g_BL <- 0.5 - 1/(1+exp(-xi))
    h_VB <- sum(h_BL + g_BL*(mq_eta-xi) - 0.5*w_BL*(Eq_sq_eta-2*xi*mq_eta+xi*xi))
    
    # Update VB log-likelihood
    Eq_loglik <- sum((y-0.5)*mq_eta) + h_VB
    
    # Update ELBO
    objval <- Eq_loglik + 0.5*Vq_logdet
    objval <- objval - lgamma(hpar$A) + lgamma(Aq_lambda)
    objval <- objval + hpar$A*log(hpar$B) - Aq_lambda*log(Bq_lambda)

    # Check convergence
    diff_abs_obj <- abs(objval-objval_old)
    diff_rel_obj <- maxreldiff(objval, objval_old)
    diff_rel_mu  <- maxreldiff(mu, mu_old)
    diff_rel_eta <- maxreldiff(mq_eta, eta_old)
    diff_rel_xi  <- maxreldiff(xi, xi_old)
    
    # Store the current state
    trace[iter,] <- c(iter, objval, Eq_lambda, diff_abs_obj, diff_rel_obj, 
                      diff_rel_mu, diff_rel_eta, diff_rel_xi)
    
    # Print the current optimization state
    if(verbose && iter%%freq==0){
      cat(gettextf(" %5d", iter),
          gettextf(" %11.5f", objval/(n+p)),
          gettextf(" %10.3e", diff_abs_obj/(n+p)),
          gettextf(" %10.3e", diff_rel_obj),
          gettextf(" %10.3e", diff_rel_mu),
          gettextf(" %10.3e", diff_rel_eta),
          gettextf(" %10.3e", diff_rel_xi),
          "\n", sep="")
    }
    
    if((diff_abs_obj<(n+p)*abstol) && (diff_rel_obj<reltol)){
      break
    }
  }
  
  # Compute the final omega
  omega = switch(solver,
                 "chol" = tcrossprod(invR),
                 "smw" = smw_chol_solve(R, X, pL2, diag(p)),
                 "sparse" = list(var=sparseinv::Takahashi_Davis(Q), cholQ=R))
  
  # Print the final optimization state
  if(verbose){
    cat(gettextf(" %5d", iter),
        gettextf(" %11.5f", objval/(n+p)),
        gettextf(" %10.3e", diff_abs_obj/(n+p)),
        gettextf(" %10.3e", diff_rel_obj),
        gettextf(" %10.3e", diff_rel_mu),
        gettextf(" %10.3e", diff_rel_eta),
        gettextf(" %10.3e", diff_rel_xi), "\n",
        paste0(rep("-", 73), collapse=""),
        "\n", sep="")
  }
  
  if (iter==maxiter)
    warning("BLVB-LOGIT: the algorithm has not reached convergence", call.=FALSE, immediate.=TRUE)
  
  
  # Set the final output
  summary$beta$mean <- mu
  summary$beta$var <- omega
  summary$beta$chol <- R
  summary$beta$sqtr <- Eq_sq_beta
  summary$beta$logdet <- Vq_logdet
  
  summary$eta$mean <- mq_eta
  summary$eta$var <- Vq_eta
  
  summary$lambda$shape <- Aq_lambda
  summary$lambda$rate <- Bq_lambda
  summary$lambda$mean <- Eq_lambda
  summary$lambda$inv <- Eq_inv_lambda
  summary$lambda$log <- Eq_log_lambda
  
  state$xi <- xi
  state$w_BL <- w_BL
  state$h_BL <- h_BL
  state$h_VB <- h_VB
  state$elbo <- objval
  state$niter <- iter
  
  return(list(summary=summary, state=state, trace=trace[1:iter,], call=setting))
}
