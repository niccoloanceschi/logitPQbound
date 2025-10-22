
## PACKAGE IMPORT ----

library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

SHOW <- FALSE
SAVE <- FALSE

# DATAPATH <- "data/Portland"
# SAVEPATH <- "tutorial/results"
# IMGPATH  <- "img/portland_RIDGE"

DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/results/Portland"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img/portland_RIDGE", sep="/")


pch <- c(15:18)
col <- c(2:4,7)

## PORTLAND DATA ----

load(paste(DATAPATH, "PortlandData.RData", sep="/"))

## UTILITY FUNCTIONS ----

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

plot_fit_path <- function (fit_list, field, log=FALSE, main="", position="bottomright", ...) {
  # Get the lambda grid
  loglambdas <- log10(fit_list[[1]]$lambdas)
  
  # Extract the specified field from each list 
  mat <- fit_list %>% lapply(\(.) .[[field]]) %>% do.call("cbind", .)
  
  # If a "norm" field is required, compute it using the estimated coefficients
  if (field== "norm") mat <- fit_list %>% lapply(\(.) sqrt(colMeans(.[["beta"]]^2))) %>% do.call("cbind", .)
  if (field=="pnorm") mat <- fit_list %>% lapply(\(.) {
    r <- abs(D %*% .[["beta"]])
    pL2 <- colMeans(r*r)/2
    return(pL2)
  }) %>% do.call("cbind", .)
  
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


## OPTIM SET-UP ----

devtools::load_all()

# Data dimensions
n <- nrow(psi)
p <- ncol(psi)

# Model matrices
y <- .5*obs+.5
X <- psi
D <- (1/sqrt(colSums(mass))) * stiff

# Penalty parameters
lambdas <- 10^seq(-6, +1, by=0.25) # for solution path
lambda <- .00001 # for single fit

# Cross-validation Seed and n of folds
seed <- 123456
nfold <- 5

# Intercept penalty
eps <- 1e-8
intercept <- FALSE

# EDF inflation factor (for GCV computation only)
gamma <- 1.

# PQ update
phi <- 0.9

# Convergence tolerance 
objtol <- 1e-7
reltol <- 1e-2
abstol <- 1e-3
etatol <- 1.
maxiter <- 1000

# Progress output
verbose <- TRUE
freq <- 10

# Initial values
beta0 <- rnorm(p, mean=0, sd=sqrt(1/p))

## SINGLE FIT ----

### BL fit ----
exetime_1run_BL <- proc.time()
fit_1run_BL <- fit_logit_spridge(y, X, D, type='BL', beta_start=beta0, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol,
                                 verbose=verbose, 
                                 freq=freq)
exetime_1run_BL <- (proc.time() - exetime_1run_BL)[3]

### PG fit ----
exetime_1run_PG <- proc.time()
fit_1run_PG <- fit_logit_spridge(y, X, D, type='PG', beta_start=beta0, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol, 
                                 verbose=verbose, freq=freq)
exetime_1run_PG <- (proc.time() - exetime_1run_PG)[3]

### PQ fit ----
exetime_1run_PQ <- proc.time()
fit_1run_PQ <- fit_logit_spridge(y, X, D, type='PQ', beta_start=beta0, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol, 
                                 verbose=verbose, freq=freq)
exetime_1run_PQ <- (proc.time() - exetime_1run_PQ)[3]

### NR fit ----
exetime_1run_NR <- proc.time()
fit_1run_NR <- fit_logit_spridge(y, X, D, type='NR', beta_start=beta0, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol, 
                                 verbose=verbose, freq=freq)
exetime_1run_NR <- (proc.time() - exetime_1run_NR)[3]

### Summary ----
fit_1run_BL$exetime <- exetime_1run_BL
fit_1run_PG$exetime <- exetime_1run_PG
fit_1run_PQ$exetime <- exetime_1run_PQ
fit_1run_NR$exetime <- exetime_1run_NR
fit_1run_list <- list(BL=fit_1run_BL, PG=fit_1run_PG, 
                      PQ=fit_1run_PQ, NR=fit_1run_NR)

df_1run_summary <- data.frame(
  method = c("BL", "PG", "PQ", "NR"),
  niter = sapply(fit_1run_list, \(.) .$niter),
  exetime = sapply(fit_1run_list, \(.) .$exetime),
  loglik = sapply(fit_1run_list, \(.) .$loglik),
  row.names = c(1:4))

print(df_1run_summary)

if (SAVE) {
  filename <- "portland_ridge_1run_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}


if (SAVE) {
  filename <- "portland_ridge_1run_loglik.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  pch <- c(15:18); col <- c(2:4,7)
  # Log-likelihood
  niter <- max(sapply(fit_1run_list, \(.) .$niter))
  objval <- matrix(NA, nrow=niter, ncol=4)
  objval[1:fit_1run_BL$niter,1] <- fit_1run_BL$trace$objval
  objval[1:fit_1run_PG$niter,2] <- fit_1run_PG$trace$objval
  objval[1:fit_1run_PQ$niter,3] <- fit_1run_PQ$trace$objval
  objval[1:fit_1run_NR$niter,4] <- fit_1run_NR$trace$objval
  colnames(objval) <- c("BL", "PG", "PQ", "NR")
  matplot(objval[2:min(30, niter),], type="b", lty=1, col=col, pch=pch, xlab="", ylab="")
  title(xlab="Iteration", ylab="Log-Likelihood", main="Penalized Log-Likelihood")
  legend("bottomright", col=col, pch=pch, legend=c("BL", "PG", "PQ", "NR"))
  # Iterations
  with(df_1run_summary, {
    barplt <- barplot(niter, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iteration", main="Number of Iterations")
  })
  # Exetime
  with(df_1run_summary, {
    barplt <- barplot(exetime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_ridge_1run_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(2,2))
    pred_1run_BL <- as.vector(1 / (1 + exp(- X %*% fit_1run_BL$beta)))
    pred_1run_PG <- as.vector(1 / (1 + exp(- X %*% fit_1run_PG$beta)))
    pred_1run_NR <- as.vector(1 / (1 + exp(- X %*% fit_1run_NR$beta)))
    pred_1run_PQ <- as.vector(1 / (1 + exp(- X %*% fit_1run_PQ$beta)))
    
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_BL, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="BL"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_PG, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PG"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_NR, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="NR"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_PQ, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PQ"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}


## SOLUTION PATH ----

### BL fit ----
exetime_path_BL <- proc.time()
fit_path_BL <- fit_logit_spridge_path(y, X, D, type='BL', beta_start=beta0, 
                                      lambda=lambdas, gamma=gamma, phi=phi, 
                                      maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                      etatol=etatol, verbose=verbose, freq=freq)
exetime_path_BL <- (proc.time() - exetime_path_BL)[3]

### PG fit ----
exetime_path_PG <- proc.time()
fit_path_PG <- fit_logit_spridge_path(y, X, D, type='PG', beta_start=beta0, 
                                      lambda=lambdas, gamma=gamma, phi=phi, 
                                      maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                      etatol=etatol, verbose=verbose, freq=freq)
exetime_path_PG <- (proc.time() - exetime_path_PG)[3]

### PQ fit ----
exetime_path_PQ <- proc.time()
fit_path_PQ <- fit_logit_spridge_path(y, X, D, type='PQ', beta_start=beta0, 
                                      lambda=lambdas, gamma=gamma, phi=phi, 
                                      maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                      etatol=etatol, verbose=verbose, freq=freq)
exetime_path_PQ <- (proc.time() - exetime_path_PQ)[3]

### NR fit ----
exetime_path_NR <- proc.time()
fit_path_NR <- fit_logit_spridge_path(y, X, D, type='NR', beta_start=beta0, 
                                      lambda=lambdas, gamma=gamma, phi=phi, 
                                      maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                      etatol=etatol, verbose=verbose, freq=freq)
exetime_path_NR <- (proc.time() - exetime_path_NR)[3]

### Summary ----

cat("BL:", exetime_path_BL, "\n",
    "PG:", exetime_path_PG, "\n",
    "PQ:", exetime_path_PQ, "\n",
    "NR:", exetime_path_NR, "\n")

fit_path_BL$tottime <- exetime_path_BL
fit_path_PG$tottime <- exetime_path_PG
fit_path_PQ$tottime <- exetime_path_PQ
fit_path_NR$tottime <- exetime_path_NR
fit_path_list <- list(BL=fit_path_BL, PG=fit_path_PG, 
                      PQ=fit_path_PQ, NR=fit_path_NR)

df_path_summary <- data.frame(
  method = c("BL", "PG", "PQ", "NR"),
  niter = sapply(fit_path_list, \(.) mean(.$niter)),
  exetime = sapply(fit_path_list, \(.) .$tottime),
  loglik = sapply(fit_path_list, \(.) mean(.$loglik)),
  row.names = c(1:4))

print(df_path_summary)

if (SAVE) {
  filename <- "portland_ridge_path_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_summary, file=filepath, row.names=FALSE)
}

df_path_extended <- data.frame(
  method = rep(c("BL", "PG", "PQ", "NR"), each=length(lambdas)),
  lambda = c(sapply(fit_path_list, \(.) .$lambda)),
  niter = c(sapply(fit_path_list, \(.) .$niter)),
  exetime = c(sapply(fit_path_list, \(.) .$exetime)),
  loglik = c(sapply(fit_path_list, \(.) .$loglik)))

print(df_path_extended)

if (SAVE) {
  filename <- "portland_ridge_path_extended.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_extended, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- "portland_ridge_path_timegain.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(df_path_summary, {
    par(mfrow=c(1,3))
    # Total number of iterations
    barplt <- barplot(29*niter, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, 29*niter-0.05*max(29*niter), round(29*niter, 2))
    title(ylab="Iterations", main="Total number of iterations")
    # Total execution time
    barplt <- barplot(exetime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Total execution Time")
    # Relative time gain
    reltime <- 100*(1-exetime[3]/exetime)
    barplt <- barplot(reltime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, reltime-0.1*sign(reltime)*max(reltime), paste(floor(reltime), "%"))
    title(ylab="Time Gain", main="Time Gain")
    par(mfrow=c(1,1))
  })
  dev.off()
}

### Solution path ----

pch <- c(15:18)
col <- c(2:4,7)

plot_fit_path(fit_path_list, field="dev", main="Deviance", position="bottomright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="logdet", main="Hessian Log-Determinant", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="reml", main="Restricted Log-Likelihood", position="bottomright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="norm", main="Parameter Norm", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="pnorm", main="Penalty Seminorm", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="niter", main="Log Number of Iterations", position="topright", log=TRUE, pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", pch=pch, col=col, lty=1)
plot_fit_path(fit_path_list, field="exetime", main="Log Execution Time", position="topright", log=TRUE, pch=pch, col=col, lty=1)


if (SAVE) {
  filename <- "portland_ridge_path_exetime.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", log=FALSE, pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", log=FALSE, pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="niter", main="Log Number of Iterations", position="topright", log=TRUE, pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Log Execution Time", position="topright", log=TRUE, pch=pch, col=col, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

if (SAVE) {
  filename <- "portland_ridge_path_loglik.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="reml", main="Hessian Log-Determinant", position="bottomright", pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="norm", main="Coefficient Norm", position="topright", pch=pch, col=col, lty=1)
  plot_fit_path(fit_path_list, field="pnorm", main="Penalty Seminorm", position="topright", pch=pch, col=col, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

### Spatial maps ----

if (SAVE) {
  filename <- "portland_ridge_path_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(2,2))
    k <- sapply(fit_path_list, \(.) which.max(.$reml))
    pred_path_BL <- c(plogis(as.vector(X %*% fit_path_BL$beta[,k["BL"]])))
    pred_path_PG <- c(plogis(as.vector(X %*% fit_path_PG$beta[,k["PG"]])))
    pred_path_NR <- c(plogis(as.vector(X %*% fit_path_NR$beta[,k["NR"]])))
    pred_path_PQ <- c(plogis(as.vector(X %*% fit_path_PQ$beta[,k["PQ"]])))
    
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_BL, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="BL"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_PG, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PG"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_NR, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="NR"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_PQ, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PQ"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

## CROSS-VALIDATION ----

### BL fit ----
exetime_cv_BL <- proc.time()
fit_cv_BL <- fit_logit_spridge_cv(y, X, D, type='BL', nfold=nfold, seed=seed, 
                                  beta_start=beta0, lambda=lambdas, gamma=gamma, 
                                  phi=phi, maxiter=maxiter, abstol=objtol, 
                                  reltol=reltol, etatol=etatol, 
                                  verbose=verbose, freq=freq)
exetime_cv_BL <- (proc.time() - exetime_cv_BL)[3]

### PG fit ----
exetime_cv_PG <- proc.time()
fit_cv_PG <- fit_logit_spridge_cv(y, X, D, type='PG', nfold=nfold, seed=seed, 
                                  beta_start=beta0, lambda=lambdas, gamma=gamma, 
                                  phi=phi, maxiter=maxiter, abstol=objtol, 
                                  reltol=reltol, etatol=etatol, 
                                  verbose=verbose, freq=freq)
exetime_cv_PG <- (proc.time() - exetime_cv_PG)[3]

### PQ fit ----
exetime_cv_PQ <- proc.time()
fit_cv_PQ <- fit_logit_spridge_cv(y, X, D, type='PQ', nfold=nfold, seed=seed, 
                                  beta_start=beta0, lambda=lambdas, gamma=gamma, 
                                  phi=phi, maxiter=maxiter, abstol=objtol, 
                                  reltol=reltol, etatol=etatol, 
                                  verbose=verbose, freq=freq)
exetime_cv_PQ <- (proc.time() - exetime_cv_PQ)[3]

### NR fit ----
exetime_cv_NR <- proc.time()
fit_cv_NR <- fit_logit_spridge_cv(y, X, D, type='NR', nfold=nfold, seed=seed, 
                                  beta_start=beta0, lambda=lambdas, gamma=gamma, 
                                  phi=phi, maxiter=maxiter, abstol=objtol, 
                                  reltol=reltol, etatol=etatol, 
                                  verbose=verbose, freq=freq)
exetime_cv_NR <- (proc.time() - exetime_cv_NR)[3]

### Summary ----

cat("BL:", exetime_cv_BL, "\n",
    "PG:", exetime_cv_PG, "\n",
    "PQ:", exetime_cv_PQ, "\n",
    "NR:", exetime_cv_NR, "\n")

fit_cv_BL$tottime <- exetime_cv_BL
fit_cv_PG$tottime <- exetime_cv_PG
fit_cv_PQ$tottime <- exetime_cv_PQ
fit_cv_NR$tottime <- exetime_cv_NR
fit_cv_list <- list(BL=fit_cv_BL, PG=fit_cv_PG, 
                    PQ=fit_cv_PQ, NR=fit_cv_NR)

df_cv_extended <- rbind(cbind(method="BL", fit_cv_BL$cv), 
                        cbind(method="PG", fit_cv_PG$cv), 
                        cbind(method="PQ", fit_cv_PQ$cv), 
                        cbind(method="NR", fit_cv_NR$cv))

print(df_cv_extended)

if (SAVE) {
  filename <- "portland_ridge_cv_extended.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_cv_extended, file=filepath, row.names=FALSE)
}

df_cv_summary <- data.frame(
  method = c("BL", "PG", "PQ", "NR"),
  niter = sapply(fit_cv_list, \(.) .$niter),
  exetime = sapply(fit_cv_list, \(.) .$tottime),
  loglik = sapply(fit_cv_list, \(.) .$loglik)
)

print(df_cv_summary)

if (SAVE) {
  filename <- "portland_ridge_cv_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_cv_summary, file=filepath, row.names=FALSE)
}


if (SAVE) {
  filename <- "portland_ridge_cv_timegain.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(df_cv_summary, {
    par(mfrow=c(1,2))
    barplt <- barplot(exetime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
    reltime <- 100*(1-exetime[3]/exetime)
    barplt <- barplot(reltime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, reltime-0.05*sign(reltime)*max(reltime), paste(floor(reltime), "%"))
    title(ylab="Time Gain", main="Time Gain")
    par(mfrow=c(1,1))
  })
  dev.off()
}

### CV solution path ----

if (SAVE) {
  filename <- "portland_ridge_cv_loglik.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  for (fit_cv_tmp in fit_cv_list){
    type <- fit_cv_tmp$type
    with(fit_cv_tmp$cv, {
      mat <- do.call(cbind, lapply(1:5, \(k) loglik[fold == k]))
      avg <- apply(mat, 1, mean)
      std <- apply(mat, 1, sd) / sqrt(5)
      matplot(log10(lambdas), mat, type="o", pch=18, lty=2, col=2:6, xlab="", ylab="", main="")
      points(log10(lambdas), avg, type="o", pch=19, col=1, lty=1, cex=1.5)
      abline(v=log10(lambdas)[which.max(avg)], col=8, lty=2)
      title(xlab=expression(log[10](lambda)), ylab="Log-Likelihood", main=type)
      legend("bottomright", pch=c(19,rep(18,5)), lty=c(1,rep(2,5)), col=c(1:6), 
             legend=c("Average", paste("Fold", 1:5)))
    })
  }
  par(mfrow=c(1,1))
  dev.off()
}



### Spatial maps ----

if (SAVE) {
  filename <- "portland_ridge_cv_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(2,2))
    pred_cv_BL <- as.vector(1 / (1 + exp(- X %*% fit_cv_BL$beta)))
    pred_cv_PG <- as.vector(1 / (1 + exp(- X %*% fit_cv_PG$beta)))
    pred_cv_PQ <- as.vector(1 / (1 + exp(- X %*% fit_cv_PQ$beta)))
    pred_cv_NR <- as.vector(1 / (1 + exp(- X %*% fit_cv_NR$beta)))
    
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pred_cv_BL, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="BL"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_cv_PG, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PG"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_cv_PQ, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PQ"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_cv_NR, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="NR"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}
