
#' @title Update beta via PQ bound optimization for ridge logistic regression
#' @keywords internal
update_enet_coord_PQ <- function(y, X, eta, w, nu, s, beta, pL1, pL2, intercept, niter=10){
  # Get the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  z <- (y-0.5-nu*s) / w
  r <- z - eta
  # Update the regression parameters via cyclical coordinate ascent
  beta <- update_enet_coord(r, X, eta, w, beta, pL1, pL2, intercept, niter)
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via PG bound optimization for ridge logistic regression
#' @keywords internal
update_enet_coord_PG <- function(y, X, eta, w, beta, pL1, pL2, intercept, niter=10){
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  z <- (y-0.5)/w
  r <- z - eta
  # Update the regression parameters via cyclical coordinate ascent
  beta <- update_enet_coord(r, X, eta, w, beta, pL1, pL2, intercept, niter)
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via BL bound optimization for ridge logistic regression
#' @keywords internal
update_enet_coord_BL <- function(y, X, eta, w, beta, pL1, pL2, intercept, niter=10){
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the sufficient statistics
  pr <- 1/(1+exp(-eta))
  w <- rep(0.25, n)
  z <- eta + (y-pr)/w
  r <- z - eta
  # Update the regression parameters via cyclical coordinate ascent
  beta <- update_enet_coord(r, X, eta, w, beta, pL1, pL2, intercept, niter)
  # Return the estimated beta
  return (beta)
}

#' @title Update beta via NR quadratic approximation for ridge logistic regression
#' @keywords internal
update_enet_coord_NR <- function(y, X, eta, w, beta, pL1, pL2, intercept, niter=10){
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Compute the probabilities, pseudo-data, and pseudo-residuals
  pr <- 1/(1+exp(-eta))
  z <- eta + (y-pr)/w
  r <- z - eta
  # Update the regression parameters via cyclical coordinate ascent
  beta <- update_enet_coord(r, X, eta, w, beta, pL1, pL2, intercept, niter)
  # Return the estimated beta
  return (beta)
}

update_enet_coord <- function(r, X, eta, w, beta, pL1, pL2, intercept, niter=10) {
  # Set the model dimensions
  n <- nrow(X)
  p <- ncol(X)
  # Set the coordinates to be updated
  if ( intercept) idx <- sample(2:p, size=p-1, replace=FALSE)
  if (!intercept) idx <- sample(1:p, size=p,   replace=FALSE)
  # Solve the linear system via iterative coordinate updates
  for (iter in 1:niter) {
    # Cycle over all the coordinates
    for (j in idx) {
      beta_j <- beta[j]
      xwx_j <- sum(w*X[,j]^2)
      xwr_j <- sum(X[,j]*w*r) + xwx_j*beta_j
      beta[j] <- soft_threshold(xwr_j, pL1[j]) / (xwx_j + pL2[j])
      delta_j <- as.vector(X[,j] * (beta[j] - beta_j))
      eta[] <- as.vector(eta + delta_j)
      r[] <- as.vector(r - delta_j)
    }
    # Update the intercept
    if (intercept) {
      delta_0 <- sum(w*r)/sum(w)
      beta[1] <- beta[1] + delta_0
      eta[] <- as.vector(eta + delta_0)
      r[] <- as.vector(r - delta_0)
    }
  }
  # Return the estimated beta
  return (beta)
}

