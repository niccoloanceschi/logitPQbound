

panel_cor <- function(x, y, ...){
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=2)[1]
  cex <- 0.3/strwidth(txt)
  text(0.5, 0.5, txt, cex=cex*r, col=2)
  text(0.5, 0.65, "Corr", cex=cex*r/2, col=1)
}

panel_linear <- function(x, y, ...) {
  # usr <- par("usr")
  # on.exit(par(usr))
  points(x, y, pch=20, col=1)
  abline(a=0, b=1, col=2)
}

plot_fit_path <- function (fit_list, field, D=NULL, log=FALSE, main="", position="bottomright", ...) {
  
  # Check the input field
  fields <- c("dev", "loglik", "logdet", "reml", "norm", "pnorm", "niter", "exetime")
  field  <- match.arg(field, fields)
  
  # Get the lambda grid
  loglambdas <- log10(fit_list[[1]]$lambdas)
  
  # Extract the specified field from each list 
  mat <- fit_list %>% lapply(\(.) .[[field]]) %>% do.call("cbind", .)
  
  # Set the penalty matrix
  if (is.null(D)) D <- diag(p)
  
  # If a "norm" field is required, compute it using the estimated coefficients
  if (field== "norm") {
    mat <- fit_list %>% lapply(\(.) {
      colMeans(.[["beta"]]^2)
    }) %>% do.call("cbind", .)
  }
  
  if (field=="pnorm") {
    mat <- fit_list %>% lapply(\(.) {
      alpha <- .[["alpha"]]
      r <- abs(D %*% .[["beta"]])
      pL1 <- colMeans(r)
      pL2 <- colMeans(r*r)/2
      return(alpha*pL1+(1-alpha)*pL2)
    }) %>% do.call("cbind", .)
  }
  
  # If specified, transform the results in the logarithmic scale
  if (log) mat <- log10(mat)
  
  # Plot the results
  matplot(loglambdas, mat, type="o", xlab="", ylab="", main="", ...)
  title(xlab=expression(log[10](lambda)), ylab=field, main=main)
  legend(position, legend=names(fit_list), ...)
  
  # If required, plot vertical lines identifying the REML/GCV/AIC/BIC solution
  if (field=="reml") {
    abline(v=loglambdas[apply(mat, 2, which.max)], lty=2, col=8)
  }
  if (field %in% c("gcv", "aic", "bic")) {
    abline(v=loglambdas[apply(mat, 2, which.min)], lty=2, col=8)
  }
}


plot_mesh <- function (mesh, incol="grey70", bndcol="grey30", ...) {
  # Set the x-y limits
  xrng <- range(mesh$nodes[, 1])
  yrng <- range(mesh$nodes[, 2])
  
  # Plot the mesh nodes
  plot(mesh$nodes, cex = 0.001, col=col, xlim=xrng, ylim=yrng, 
       xlab="", ylab="", xaxt="n", yaxt="n", bty="n", ...)
  
  # Plot the mesh internal segments
  x_in_start <- mesh$nodes[mesh$edges[, 1], 1]
  y_in_start <- mesh$nodes[mesh$edges[, 2], 1]
  x_in_end <- mesh$nodes[mesh$edges[, 1], 2]
  y_in_end <- mesh$nodes[mesh$edges[, 2], 2]
  segments(x_in_start, x_in_end, y_in_start, y_in_end, col=incol)
  
  # Plot the mesh boundary segments
  x_bnd_start <- mesh$nodes[mesh$segments[, 1], 1]
  y_bnd_start <- mesh$nodes[mesh$segments[, 2], 1]
  x_bnd_end <- mesh$nodes[mesh$segments[, 1], 2]
  y_bnd_end <- mesh$nodes[mesh$segments[, 2], 2]
  segments(x_bnd_start, x_bnd_end, y_bnd_start, y_bnd_end, col=bndcol, lwd=1.5)
}

plot_field <- function(coefs, basis, ngrid, ...) {
  # Set the x-y limits
  xrng <- range(basis$mesh$nodes[,1])
  yrng <- range(basis$mesh$nodes[,2])
  
  # Set the x-y grid
  xs <- seq(from=xrng[1], to=xrng[2], length=ngrid)
  ys <- seq(from=yrng[1], to=yrng[2], length=ngrid)
  
  # Set the x-y lattice
  xx <- rep(xs, ngrid)
  yy <- rep(ys, rep(ngrid, ngrid))
  
  # Set the FEM object
  fem <- fdaPDE::FEM(coefs, basis)
  
  # Compute the FEM surface
  field <- fdaPDE::eval.FEM(fem, cbind(xx, yy))
  field <- matrix(field, nrow=ngrid, ncol=ngrid)
  
  # Plot the field
  plot3D::image2D(x=xs, y=ys, z=field, colvar=field, xlab="", ylab="", ...)
}

## plot_fit_path <- function (fit_list, field, log=FALSE, main="", position="bottomright", ...) {
##   # Get the lambda grid
##   loglambdas <- log10(fit_list[[1]]$lambdas)
##   
##   # Extract the specified field from each list 
##   mat <- fit_list %>% lapply(\(.) .[[field]]) %>% do.call("cbind", .)
##   
##   # If a "norm" field is required, compute it using the estimated coefficients
##   if (field== "norm") mat <- fit_list %>% lapply(\(.) sqrt(colMeans(.[["beta"]]^2))) %>% do.call("cbind", .)
##   if (field=="pnorm") mat <- fit_list %>% lapply(\(.) {
##     r <- abs(D %*% .[["beta"]])
##     pL2 <- colMeans(r*r)/2
##     return(pL2)
##   }) %>% do.call("cbind", .)
##   
##   # If specified, transform the results in the logarithmic scale
##   if (log) mat <- log10(mat)
##   
##   # Plot the results
##   matplot(loglambdas, mat, type="o", xlab="", ylab="", main="", ...)
##   title(xlab=expression(log[10](lambda)), ylab=field, main=main)
##   legend(position, legend=names(fit_list), ...)
##   
##   # If required, plot vertical lines identifying the REML/GCV/AIC/BIC solution
##   if (field=="reml") {
##     abline(v=loglambdas[apply(mat, 2, which.max)], lty=2, col=8)
##   }
##   if (field %in% c("gcv", "aic", "bic")) {
##     abline(v=loglambdas[apply(mat, 2, which.min)], lty=2, col=8)
##   }
## }
