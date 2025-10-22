
## PACKAGE IMPORT ----

library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

SHOW <- FALSE
SAVE <- FALSE
DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/results"
IMGPATH  <- "img/portland_MFVB"

## PORTLAND DATA ----

load(paste(DATAPATH, "PortlandData.RData", sep="/"))

## UTILITY FUNCTIONS ----

plot_mesh <- function(mesh, incol="grey70", bndcol="grey30", ...) {
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

accuracy <- function(sim, mu, sigma2) {
  n <- length(sim)
  K <- 201
  idx <- seq(from=floor(n/2), to=n, by=5)
  kde <- density(sim[idx])
  lower <- min(kde$x)
  upper <- max(kde$x)
  x <- seq(from=lower, to=upper, length=K)
  f1 <- approxfun(kde$x, kde$y, rule=2)(x)
  f2 <- dnorm(x, mean=mu, sd=sqrt(sigma2))
  dx <- diff(x)
  df <- abs(f1 - f2)
  err <- 0.5*sum((df[-1]+df[-K])*dx)
  acc <- 1 - 0.5*err
  return(acc)
}

## VB-PDE FIT ----

devtools::load_all()

n <- nrow(psi)
p <- ncol(psi)
y <- .5*obs+.5
X <- psi
D <- (1/sqrt(colSums(mass))) * stiff
lambda <- 0.00001
eps <- 1e-8
intercept <- FALSE
abstol <- 1e-10
reltol <- 1e-8
maxiter <- 1000L
verbose <- TRUE
freq <- 10L

beta0 <- rnorm(p, mean=0, sd=sqrt(1/p))
beta0 <- NULL

### ---- BL-VB ----
exetime_BL <- proc.time()
fit_BL <- fit_logit_mfvb(y=y, X=X, D=D, type="BL", beta_start=beta0, 
                         lambda=lambda, eps=eps, intercept=intercept, 
                         solver="sparse", maxiter=maxiter, abstol=abstol, 
                         reltol=reltol, verbose=verbose, freq=freq)
exetime_BL <- (proc.time() - exetime_BL)[3]

### ---- PG-VB ----
exetime_PG <- proc.time()
fit_PG <- fit_logit_mfvb(y=y, X=X, D=D, type="PG", beta_start=beta0, 
                         lambda=lambda, eps=eps, intercept=intercept, 
                         solver="sparse", maxiter=maxiter, abstol=abstol, 
                         reltol=reltol, verbose=verbose, freq=freq)
exetime_PG <- (proc.time() - exetime_PG)[3]

### ---- PQ-VB ----
exetime_PQ <- proc.time()
fit_PQ <- fit_logit_mfvb(y=y, X=X, D=D, type="PQ", beta_start=beta0, 
                         lambda=lambda, eps=eps, intercept=intercept, 
                         solver="sparse", maxiter=maxiter, abstol=abstol, 
                         reltol=reltol, verbose=verbose, freq=freq)
exetime_PQ <- (proc.time() - exetime_PQ)[3]

### ---- MCMC ----
exetime_MC <- proc.time()
fit_MC <- fit_logit_mcmc(y=y, X=X, D=D, solver="sparse",
                         lambda=lambda, eps=eps, intercept=intercept,
                         maxiter=5000L, verbose=verbose, freq=100L)
exetime_MC <- (proc.time() - exetime_MC)[3]


## RESULTS ----

### Summary measures ----

# Number of iterations
niter_BL <- fit_BL$niter
niter_PG <- fit_PG$niter 
niter_PQ <- fit_PQ$niter 
niter_MC <- 5000L

# ELBO
elbo_BL <- fit_BL$elbo
elbo_PG <- fit_PG$elbo 
elbo_PQ <- fit_PQ$elbo 
elbo_MC <- NaN

# Posterior means
mu_BL <- fit_BL$mu
mu_PG <- fit_PG$mu
mu_PQ <- fit_PQ$mu
mu_MC <- fit_MC$mbeta

# Posterior variances
var_BL <- diag(fit_BL$omega$var)
var_PG <- diag(fit_PG$omega$var)
var_PQ <- diag(fit_PQ$omega$var)
var_MC <- fit_MC$vbeta

# Posterior covariance matrices
cov_BL <- as.matrix(Matrix::solve(fit_BL$omega$cholQ))
cov_PG <- as.matrix(Matrix::solve(fit_PG$omega$cholQ))
cov_PQ <- as.matrix(Matrix::solve(fit_PQ$omega$cholQ))
cov_MC <- cov(fit_MC$beta[seq(from=2501L, to=5000L, by=5),])

# Posterior linear predictor means
eta_BL <- fit_BL$eta
eta_PG <- fit_PG$eta
eta_PQ <- fit_PQ$eta
eta_MC <- fit_MC$meta

# Posterior linear predictor variances
var_eta_BL <- fit_BL$veta
var_eta_PG <- fit_PG$veta
var_eta_PQ <- fit_PQ$veta
var_eta_MC <- fit_MC$veta

# Posterior linear predictor covariance
cov_eta_BL <- X %*% tcrossprod(cov_BL, X)
cov_eta_PG <- X %*% tcrossprod(cov_PG, X)
cov_eta_PQ <- X %*% tcrossprod(cov_PQ, X)
cov_eta_MC <- X %*% tcrossprod(cov_MC, X)

# Posterior linear predictor means
pr_BL <- 1 / (1 + exp(-eta_BL))
pr_PG <- 1 / (1 + exp(-eta_PG))
pr_PQ <- 1 / (1 + exp(-eta_PQ))
pr_MC <- 1 / (1 + exp(-eta_MC))

# Accuracy scores (basis parameters)
acc_beta_BL <- sapply(1:p, function(j) {accuracy(fit_MC$beta[,j], mu_BL[j], var_BL[j])})
acc_beta_PG <- sapply(1:p, function(j) {accuracy(fit_MC$beta[,j], mu_PG[j], var_PG[j])})
acc_beta_PQ <- sapply(1:p, function(j) {accuracy(fit_MC$beta[,j], mu_PQ[j], var_PQ[j])})
acc_beta_MC <- rep(1.0, times=p)

# Accuracy scores (linear predictor)
acc_eta_BL <- sapply(1:n, function(i) {accuracy(c(fit_MC$beta %*% X[i,]), eta_BL[i], var_eta_BL[i])})
acc_eta_PG <- sapply(1:n, function(i) {accuracy(c(fit_MC$beta %*% X[i,]), eta_PG[i], var_eta_PG[i])})
acc_eta_PQ <- sapply(1:n, function(i) {accuracy(c(fit_MC$beta %*% X[i,]), eta_PQ[i], var_eta_PQ[i])})
acc_eta_MC <- rep(1.0, times=n)

# Posterior mean error (basis parameters)
err_mean_beta_BL <- abs(mu_BL-mu_MC)
err_mean_beta_PG <- abs(mu_PG-mu_MC)
err_mean_beta_PQ <- abs(mu_PQ-mu_MC)
err_mean_beta_MC <- rep(0.0, times=p)

# Posterior log-variance error (basis parameters)
err_var_beta_BL <- abs(log(var_BL)-log(var_MC))
err_var_beta_PG <- abs(log(var_PG)-log(var_MC))
err_var_beta_PQ <- abs(log(var_PQ)-log(var_MC))
err_var_beta_MC <- rep(0.0, times=p)

# Posterior mean error (linear predictor)
err_mean_eta_BL <- abs(eta_BL-eta_MC)
err_mean_eta_PG <- abs(eta_PG-eta_MC)
err_mean_eta_PQ <- abs(eta_PQ-eta_MC)
err_mean_eta_MC <- rep(0.0, times=n)

# Posterior standard deviation error (linear predictor)
err_var_eta_BL <- abs(log(var_eta_BL)-log(var_eta_MC))
err_var_eta_PG <- abs(log(var_eta_PG)-log(var_eta_MC))
err_var_eta_PQ <- abs(log(var_eta_PQ)-log(var_eta_MC))
err_var_eta_MC <- rep(0.0, times=n)

# Average accuracy/error scores
tvd_pdf_beta  <- rowMeans(rbind(1-acc_beta_BL, 1-acc_beta_PG, 1-acc_beta_PQ, 1-acc_beta_MC))
tvd_pdf_eta   <- rowMeans(rbind(1-acc_eta_BL, 1-acc_eta_PG, 1-acc_eta_PQ, 1-acc_eta_MC))
mae_mean_beta <- rowMeans(rbind(err_mean_beta_BL, err_mean_beta_PG, err_mean_beta_PQ, err_mean_beta_MC))
mae_var_beta  <- rowMeans(rbind(err_var_beta_BL, err_var_beta_PG, err_var_beta_PQ, err_var_beta_MC))
mae_mean_eta  <- rowMeans(rbind(err_mean_eta_BL, err_mean_eta_PG, err_mean_eta_PQ, err_mean_eta_MC))
mae_var_eta   <- rowMeans(rbind(err_var_eta_BL, err_var_eta_PG, err_var_eta_PQ, err_var_eta_MC))

# Summary table
summary <- data.frame(
  method = c("BL", "PG", "PQ", "MC"),
  niter = c(niter_BL, niter_PG, niter_PQ, niter_MC),
  exetime = c(exetime_BL, exetime_PG, exetime_PQ, exetime_MC),
  elbo = c(elbo_BL, elbo_PG, elbo_PQ, elbo_MC),
  tv_dist_beta = round(tvd_pdf_beta, 3),
  mae_mean_beta = round(mae_mean_beta, 3),
  mae_var_beta = round(mae_var_beta, 3),
  tv_dist_eta = round(tvd_pdf_eta, 3),
  mae_mean_eta = round(mae_mean_eta, 3),
  mae_var_eta = round(mae_var_eta, 3),
  row.names = 1:4
)

print(summary)

if (SAVE) {
  filename <- "portland_mfvb_summary.csv"
  filepath <- paste(SAVEPATH, filename, sep="/")
  write.csv2(summary, file=filepath, row.names=FALSE)
}

### summary plots ----

pch <- c(15:18)
col <- c(2:4,7)

if (SAVE) {
  filename <- "portland_mfvb_elbo.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  # ELBO
  xlim <- c(0, max(niter_BL, niter_PG, niter_PQ))
  ylim <- range(c(fit_BL$trace$objval, fit_PG$trace$objval, fit_PQ$trace$objval))/n
  plot(0, 0, pch=19, col=3, xlim=xlim, ylim=ylim, type="o", xlab="", ylab="")
  points(fit_BL$trace$iter, fit_BL$trace$objval/n, pch=pch[1], col=col[1], type="o")
  points(fit_PG$trace$iter, fit_PG$trace$objval/n, pch=pch[2], col=col[2], type="o")
  points(fit_PQ$trace$iter, fit_PQ$trace$objval/n, pch=pch[3], col=col[3], type="o")
  title(xlab="Iterations", ylab="ELBO", main="Evidence lower bound")
  legend("bottomright", pch=c(15,19,17), col=2:4, legend=c("BL-VB", "PG-VB", "PQ-VB"))
  # n-iter
  with(summary[1:3,], {
    barplt <- barplot(niter, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iteration", main="Number of Iterations")
  })
  # exe-time
  with(summary[1:3,], {
    barplt <- barplot(exetime, names.arg=method, col=col, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_error_beta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(1,3))
  # TV distance
  df <- cbind("BL-VB"=1-acc_beta_BL, "PG-VB"=1-acc_beta_PG, "PQ-VB"=1-acc_beta_PQ)
  boxplot(df, col=c(2,3,4), ylim=c(0,1), main=expression(d[TV](q[VB], p[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - mu
  df <- cbind("BL-VB"=err_mean_beta_BL, "PG-VB"=err_mean_beta_PG, "PQ-VB"=err_mean_beta_PQ)
  boxplot(df, col=c(2,3,4), main=expression(MAE(mu[VB], mu[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - log(sigma)
  df <- cbind("BL-VB"=err_var_beta_BL, "PG-VB"=err_var_beta_PG, "PQ-VB"=err_var_beta_PQ)
  boxplot(df, col=c(2,3,4), main=expression(MAE(log~sigma[VB], log~sigma[MC])))
  abline(h=.0, col=7, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_error_eta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(1,3))
  # TV distance
  boxplot(cbind("BL-VB"=1-acc_eta_BL, "PG-VB"=1-acc_eta_PG, "PQ-VB"=1-acc_eta_PQ), 
          col=c(2,3,4), ylim=c(0,1), main=expression(d[TV](q[VB], p[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - mu
  boxplot(cbind("BL-VB"=err_mean_eta_BL, "PG-VB"=err_mean_eta_PG, "PQ-VB"=err_mean_eta_PQ), 
          col=c(2,3,4), main=expression(MAE(mu[VB], mu[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - log(sigma)
  boxplot(cbind("BL-VB"=err_var_eta_BL, "PG-VB"=err_var_eta_PG, "PQ-VB"=err_var_eta_PQ), 
          col=c(2,3,4), main=expression(MAE(log~sigma[VB], log~sigma[MC])))
  abline(h=.0, col=7, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

### Spatial maps ----

if (SAVE) {
  filename <- "portland_mfvb_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pr_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_MC, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="MCMC"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}
  
if (SAVE) {
  filename <- "portland_mfvb_map_eta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(eta_BL, eta_PG, eta_PQ, eta_MC))
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=eta_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=eta_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=eta_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=eta_MC, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="MCMC"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_map_var.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(sqrt(cbind(var_eta_BL, var_eta_PG, var_eta_PQ, var_eta_MC)))
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=sqrt(var_eta_BL), xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=sqrt(var_eta_PG), xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=sqrt(var_eta_PQ), xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=sqrt(var_eta_MC), xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="MCMC"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_map_tvd.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1,3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=1-acc_eta_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=1-acc_eta_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=1-acc_eta_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_map_mae_mu.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1,3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(err_mean_eta_BL, err_mean_eta_PG, err_mean_eta_PQ))
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=err_mean_eta_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=err_mean_eta_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=err_mean_eta_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_map_mae_var.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1,3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(err_var_eta_BL, err_var_eta_PG, err_var_eta_PQ))
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=err_var_eta_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=err_var_eta_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=err_var_eta_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

### Spatial fields ----

if (SAVE) {
  filename <- "portland_mfvb_field_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    plot_field(c(plogis(mu_BL)), basis, 200, clim=c(0,1), main="BL-VB")
    plot_field(c(plogis(mu_PG)), basis, 200, clim=c(0,1), main="PG-VB")
    plot_field(c(plogis(mu_PQ)), basis, 200, clim=c(0,1), main="PQ-VB")
    plot_field(c(plogis(mu_MC)), basis, 200, clim=c(0,1), main="MCMC")
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_field_eta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(mu_BL, mu_PG, mu_PQ, mu_MC))
    plot_field(mu_BL, basis, 200, clim=clim, main="BL-VB")
    plot_field(mu_PG, basis, 200, clim=clim, main="PG-VB")
    plot_field(mu_PQ, basis, 200, clim=clim, main="PQ-VB")
    plot_field(mu_MC, basis, 200, clim=clim, main="MCMC")
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_field_var.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(sqrt(cbind(var_BL, var_PG, var_PQ)))
    plot_field(sqrt(var_BL), basis, 200, clim=clim, main="BL-VB")
    plot_field(sqrt(var_PG), basis, 200, clim=clim, main="PG-VB")
    plot_field(sqrt(var_PQ), basis, 200, clim=clim, main="PQ-VB")
    plot_field(sqrt(var_MC), basis, 200, clim=clim, main="MCMC")
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_field_tvd.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1,3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(1-acc_beta_BL, 1-acc_beta_PG, 1-acc_beta_PQ, 1-acc_beta_MC))
    plot_field(1-acc_beta_BL, basis, 200, clim=clim, main="BL-VB")
    plot_field(1-acc_beta_PG, basis, 200, clim=clim, main="PG-VB")
    plot_field(1-acc_beta_PQ, basis, 200, clim=clim, main="PQ-VB")
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_field_mae_mu.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1, 3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(err_mean_beta_BL, err_mean_beta_PG, err_mean_beta_PQ, err_mean_beta_MC))
    plot_field(err_mean_beta_BL, basis, 200, clim=clim, main="BL-VB")
    plot_field(err_mean_beta_PG, basis, 200, clim=clim, main="PG-VB")
    plot_field(err_mean_beta_PQ, basis, 200, clim=clim, main="PQ-VB")
    par(mfrow=c(1,1))
  })
  dev.off()
}

if (SAVE) {
  filename <- "portland_mfvb_field_mae_var.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(1,3))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
    clim <- range(cbind(err_var_beta_BL, err_var_beta_PG, err_var_beta_PQ, err_var_beta_MC))
    plot_field(err_var_beta_BL, basis, 200, clim=clim, main="BL-VB")
    plot_field(err_var_beta_PG, basis, 200, clim=clim, main="PG-VB")
    plot_field(err_var_beta_PQ, basis, 200, clim=clim, main="PQ-VB")
    par(mfrow=c(1,1))
  })
  dev.off()
}

## PAPER PLOT ----

ggplot_field_grid <- function(coefs, basis, ngrid, locs=NULL, ...) {
  H <- ncol(coefs)
  
  # Set the x-y limits
  xrng <- range(basis$mesh$nodes[,1])
  yrng <- range(basis$mesh$nodes[,2])
  
  # Set the x-y grid
  xs <- seq(from=xrng[1], to=xrng[2], length=ngrid)
  ys <- seq(from=yrng[1], to=yrng[2], length=ngrid)
  
  # Set the x-y lattice
  xx <- rep(xs, ngrid)
  yy <- rep(ys, rep(ngrid, ngrid))
  
  # Compute the FEM surface
  field <- matrix(NA, nrow=ngrid^2, ncol=ncol(coefs))
  for(h in 1:H){
    fem <- fdaPDE::FEM(coefs[,h], basis)
    field[,h] <- as.vector(fdaPDE::eval.FEM(fem, cbind(xx, yy)))
  }
  
  # Set the data-frames
  df_fems <- data.frame(
    g = rep(colnames(coefs), each=ngrid^2),
    x = rep(xx, times=H),
    y = rep(yy, times=H),
    z = as.vector(field))
  
  df_txt <- data.frame(
    g = colnames(coefs),
    x = rep(xrng[2]-.150*diff(xrng), H),
    y = rep(yrng[2]-.025*diff(yrng), H),
    t = paste0("TV = ", round(colMeans(field, na.rm=TRUE), 4)))
  
  df_locs <- as.data.frame(locs)
  df_locs$z <- factor(df_locs$z, labels=c("safe","risky"))
  
  # Create the data-frame
  plt <- ggplot() + theme_bw() +
    geom_raster(data=df_fems, map=aes(x=x, y=y, fill=z)) + 
    geom_label(data=df_txt, map=aes(x=x, y=y, label=t)) + 
    facet_grid(cols=vars(g)) +
    scale_shape_manual(values=c(1,19)) +
    scale_fill_viridis_c(limits=c(0,1), na.value="transparent") +
    labs(x="Longitude", y="Latitude", shape="Response", fill="TV")
  
  # Plot the field
  return(plt)
}

if (SAVE) {
  filename <- "portland_mfvb_field_tvd_gg.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.25
  tvd <- 1-cbind(BL=acc_beta_BL, PG=acc_beta_PG, PQ=acc_beta_PQ)
  dat <- data.frame(x=locs[,1], y=locs[,2], z=obs)
  plt <- ggplot_field_grid(coefs=tvd, basis=basis, ngrid=200, locs=dat)
  plt <- plt + theme(axis.title=element_blank(), axis.ticks=element_blank(), 
                     axis.text=element_blank(), strip.text=element_text(size=15),
                     legend.key.height=unit(.125, units="npc"),
                     legend.key.width=unit(.0125, units="npc"))
  ggplot2::ggsave(filename=filename, plot=plt, path=IMGPATH, 
                  width=zoom*width, height=zoom*height)
}

## END OF FILE ----
