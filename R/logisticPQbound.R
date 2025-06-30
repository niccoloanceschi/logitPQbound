
#' accurate inverse-matrix vector product
#' 
#' @param cholQ (real matrix) Cholesky matrix
#' @param r (real vector) vector to be pre-multiplied with inverse
#'
#' @export
#' 
invQ_r <- function(cholQ,r){
  return(backsolve(cholQ,forwardsolve(t(cholQ),r)))
}

#' logistic penalized log likelihood under elastic net regularization
#' 
#' @import stats
#' @importFrom Rcpp evalCpp sourceCpp
#' @useDynLib logitPQbound
#'
#' @param eta (real vector) linear predictors
#' @param y (binary vector) responses
#' @param beta (real vector) parameters (including intercept)
#' @param alpha (real) relative weight of L1 and L2 penalty
#' @param lambda (real) regularization parameter
#' 
#' @return logistic penalized log likelihood
#' 
#' @export
logpost_R <- function(eta,y,beta,alpha,lambda){
  logpost(eta,y,beta,alpha,lambda)
}


#' Ridge logistic regression via MM optimization
#'
#'
#' @param x  A \verb{n x p} dimensional design matrix. The intercept must be manually included.
#' @param y  A \verb{n} dimensional binary vector.
#' @param type (string) Family of surrogate functions for iterative optimization. Must be one of 'NR', 'BL', 'PG', or 'PQ'
#' @param beta_start (real vector) Initialization point. If missing, the algorithm is initialized setting at zero all the coefficients.
#' @param lambda (real) Regularization parameters for ridge regression
#' @param maxiter (integer) Maximum number of iterations. If reached, the algorithm raises an error.
#' @param tol (real) Threshold at which the algorithm stops.
#' @param eps_pen (real) Small value added to the diagonal of the Hessian matrix to avoid numerical issues.
#' @param maxiter_in (integer) Maximum number of inner iterations for PQ bound optimization (either for projected gradient or quadratic programming).
#' @param tol_in (real) Threshold at which the inner optimization of PQ bound stops (either for projected gradient or quadratic programming).
#' @param alpha_prjg (real) Step size for projected gradient descent on PQ bounds.
#' @param use_qp (boolean) Use quadratic programming inner optimization of PQ bounds.
#' @param use_nn (boolean) Use nearest neighbors search as safety check on PQ bounds optimization in later iterations.
#' 
#' @return The function returns a list with optimal coefficents and likelihood path through MM iterations
#' 
#' @export
#' 
ridge_logit <- function(x, y, type='PQ', beta_start=NULL, lambda=NULL,
                        maxiter=1e4, tol=1e-10, eps_pen=1e-10, 
                        maxiter_in=1e2, tol_in=1e-1, alpha_prjg=0.1, 
                        use_qp=F, use_nn=F){
  
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  
  if(!type%in%c('NR','BL','PG','PQ')){
    stop("Type must be one of 'NR', 'BL', 'PG', or 'PQ'")
  }
  
  p <- ncol(x)
  n <- nrow(x)
  
  loglik <- numeric(maxiter)
  
  # Auxiliary variables ----
  
  pL2 <- rep(0.,p)
  if(p>n){pL2 <- rep(eps_pen,p)}
  if(!is.null(lambda) & lambda!=0){pL2[-1] <- lambda}
  
  if(p>n){
    x_pL2 <- t(t(x)/pL2)
    xPxt  <- x%*%t(x_pL2)
    ptXy  <- crossprod(x_pL2,y-0.5)
  } else {
    xy <- crossprod(x,y-0.5)
  }
  
  # Bound initialization ----
  
  t0=1
  
  if(is.null(beta_start)){beta_start = rep(0,p)}
  beta <- beta_start
  
  eta  <- c(x%*%beta)
  
  if (type=='PQ') {
    eta_stable = F
    if(all(beta==0)){
      loglik[t0] <- logpost_R(eta,y,beta,0.,lambda)
      t0 <- t0+1
      if(p>n){
        if(use_qp){qp_pars = qpmadr::qpmadParameters(maxIter=maxiter_in,tol=tol_in)}
        cholQ <- chol(diag(4,n,n)+xPxt)
        d0_qp <- xPxt%*%(y-0.5)
        beta  <- crossprod(x_pL2,y-0.5) - crossprod(x_pL2,invQ_r(cholQ,d0_qp))
        eta <- c(x%*%beta)
      } else {
        cholQ <- chol(diag(pL2,p,p)+crossprod(x,0.25*x))
        beta  <- invQ_r(cholQ,xy)
        eta <- c(x%*%beta)
      }
    }
    coeff_pq <- coeff_PQ(eta)
    w  <- coeff_pq[1:n]
    nu <- coeff_pq[(n+1):(2*n)]
  } else if (type=='PG'){
    w  <- c(w_PG(eta))
  } else if (type=='NR'){      
    prob <- 1/(1+exp(-eta))
    w  <- prob*(1-prob)
  } else { # if (type=='BL')
    prob <- 1/(1+exp(-eta))
    w  <- 0.25
    if(p>n){
      cholQ <- chol(xPxt + diag(1/w,n,n))
    } else {
      cholQ <- chol(diag(pL2,p,p)+crossprod(x,w*x))
    }
  } 
  
  loglik[t0] <- logpost_R(eta,y,beta,0.,lambda)
  t0 <- t0+1
  
  # Iterative procedure ----
  
  for(t in t0:maxiter){
    
    ## Update beta ----
    
    if(p>n){
      if (type=='PQ'){
        cholQ   <- chol(xPxt + diag(1/w,n,n))
        xPxt_nu <- t(nu*xPxt)
        u_qp = sign(eta)
        if(!eta_stable){
          if(use_qp){ ## PQ optim via quadratic programming ------------------ #
            C_qp <- nu*xPxt_nu - crossprod(xPxt_nu,invQ_r(cholQ,xPxt_nu))
            d_qp <- nu*d0_qp - crossprod(xPxt_nu,invQ_r(cholQ,d0_qp))
            u_qp = qpmadr::solveqp(C_qp,h=-d_qp,lb=-1,ub=1,pars=qp_pars)$solution
          } else {    ## PQ optim via projected descent ---------------------- #
            u_qp = prjgrad_large_p(u_qp, (y-0.5)/nu, nu/w, xPxt_nu, cholQ, alpha_prjg, tol_in, maxiter_in)
            u_qp = sign(as.vector(u_qp))
          }
        } else if (use_nn) { ## PQ optim via nearest-neighbors vertex search - #
          u_vertex = u_qp%x%t(rep(1,n))
          diag(u_vertex) = -diag(u_vertex)
          u_vertex = cbind(u_qp,u_vertex)
          u_diff = u_vertex - ((y-0.5)/nu)%x%t(rep(1,n+1))
          f_vertex = 0.5*colSums(u_diff*crossprod(xPxt_nu,invQ_r(cholQ,(nu/w)*u_diff)))
          u_qp = u_vertex[,which.min(f_vertex)[1]]
        }
        beta <- crossprod(x_pL2,y-0.5-nu*u_qp) - crossprod(x_pL2,invQ_r(cholQ,xPxt%*%(y-0.5-nu*u_qp)))
      } else if (type=='PG'){
        cholQ <- chol(xPxt+diag(1/w,n,n))
        beta  <- ptXy - crossprod(x_pL2,invQ_r(cholQ,xPxt%*%(y-0.5)))
      } else {
        if (type=='NR'){cholQ <- chol(xPxt+diag(1/w,n,n))}
        beta <- crossprod(x_pL2,y-prob+w*eta) - crossprod(x_pL2,invQ_r(cholQ,xPxt%*%(y-prob+w*eta)))
      }
    } else {
      if (type=='PQ'){
        cholQ <- chol(diag(pL2,p,p)+crossprod(x,w*x))
        nu_x  <- nu*x
        u_qp = sign(eta)
        if(!eta_stable){
          if(F){   ## PQ optim via quadratic programming --------------------- #
            C_qp <- nu_x%*%invQ_r(cholQ,t(nu_x)) + diag(eps_pen,n,n)
            d_qp <- nu_x%*%invQ_r(cholQ,xy)
            u_qp = qpmadr::solveqp(C_qp,h=-d_qp,lb=-1,ub=1,pars=qp_pars)$solution
          } else { ## PQ optim via projected descent ------------------------- #
            u_qp = prjgrad_large_n(u_qp, (y-0.5)/nu, nu_x, cholQ, alpha_prjg, tol_in, maxiter_in)
            u_qp = sign(as.vector(u_qp))
          }
        } else if (use_nn) { ## PQ optim via nearest-neighbors vertex search - #
          XtNut = t(u_qp*nu_x)
          QXtNut = invQ_r(cholQ,XtNut)
          uNXQXtNut = colSums(XtNut*QXtNut)
          f_vertex = c(0,2*uNXQXtNut-2*as.vector(t(rowSums(XtNut)-xy)%*%QXtNut))
          idx_min = which.min(f_vertex)[1]
          u_qp[idx_min-1] = -u_qp[idx_min-1]
        }
        beta <- invQ_r(cholQ,xy-crossprod(x,u_qp*nu))
      } else if (type=='PG'){
        cholQ <- chol(diag(pL2,p,p)+crossprod(x,w*x))
        beta  <- invQ_r(cholQ,xy)
      } else {
        if (type=='NR'){cholQ <- chol(diag(pL2,p,p)+crossprod(x,w*x))}
        beta  <- beta + invQ_r(cholQ,crossprod(x,y-prob)-pL2*beta)
      }
    }
    
    # Update bound ---- 
    
    etaOLD <- eta
    eta <- c(x%*%beta)
    
    if (type=='PQ'){         
      coeff_pq <- coeff_PQ(eta)
      w  <- coeff_pq[1:n]
      nu <- coeff_pq[(n+1):(2*n)]
      eta_stable <- (max(abs(eta-etaOLD)/abs(etaOLD))<tol_in) & (min(eta*etaOLD)>0)
    } else if (type=='PG'){ 
      w  <- c(w_PG(eta))
    } else { 
      prob <- 1/(1+exp(-eta))
      if (type=='NR') {w <- prob*(1-prob)}
    }
    
    # Check convergence ----
    
    loglik[t] <- logpost_R(eta,y,beta,0.,lambda)
    
    if(abs(loglik[t] - loglik[t-1]) < tol){
      return(list(beta=beta,Convergence=cbind(Iteration=(1:t)-1, Loglikelihood=loglik[1:t])))
    }
  }
  stop("The algorithm has not reached convergence")
}




  
  
  
  
  
  
  
  
  
  
  
  