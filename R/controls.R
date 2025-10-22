
#' @title Set the control parameters of the projected gradient algorithm
#' @keywords internal
set_ctr_prjg = function(stepsize=1e-1, tol=1e-2, maxiter=100){
  list(stepsize=stepsize, tol=tol, maxiter=maxiter)
}

#' @title Set the control parameters of the dual QP algorithm
#' @keywords internal
set_ctr_dual = function(tol=1e-4, maxiter=100){
  list(tol=tol, maxiter=maxiter)
}

#' @title Set the control parameters of the ADMM algorithm
#' @keywords internal
set_ctr_admm = function(rho=1.0, tau=0.0, gamma=0.0, intercept=FALSE, 
                        spthr=0.9, smw=FALSE, precondition=FALSE,
                        objtol=1e-3, reltol=1e-3, abstol=1e-4, 
                        maxiter=100, verbose=FALSE, freq=10){
  list(rho=rho, tau=tau, gamma=gamma, intercept=intercept, 
       spthr=spthr, smw=smw, precondition=precondition,
       objtol=objtol, reltol=reltol, abstol=abstol, 
       maxiter=maxiter, verbose=verbose, freq=freq)
}
