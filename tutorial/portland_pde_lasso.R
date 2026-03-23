
## PACKAGE IMPORT ----

library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

SHOW <- TRUE
SAVE <- TRUE

DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/results/Portland"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img/portland_LASSO", sep="/")

DATALAB <- "portland"

# pch <- c(15:18)
# col <- c(2:4,7)

COLORS  <- c(2:4,7)
MARKERS <- c(15:18)

## PORTLAND DATA ----

load(paste(DATAPATH, "PortlandData.RData", sep="/"))

## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

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
alphas <- seq(.05, 1, by=0.15) # for solution path

lambda <- .0001 # for single fit
alpha <- alphas[6] # Mixing weight of L1 and L2 penalties

lalpha <- as.character(100*alpha)
lalpha <- ifelse(100*alpha>10, lalpha, paste0("0", lalpha))

# Cross-validation Seed and n of folds
seed <- 123456
nfold <- 5

# Intercept penalty
eps <- 1e-8
intercept <- FALSE

# EDF inflation factor (for GCV computation only)
gamma <- 1.

# PQ approximate update
phi <- 0.5
approx <- FALSE

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

# ADMM control parameters
ctr <- set_ctr_admm(
  rho=sqrt(p/n)*lambda, # Augmented Lagrangian penalty parameter
  gamma=0.01, # Proximal penalty parameter
  smw=FALSE, # Sherman-Morrison-Woodbury update
  precondition=FALSE, # Preconditioned ADMM
  objtol=.1*objtol, # Tolerance for the objective function
  reltol=reltol, # Tolerance for the relative change of the primal-dual residuals
  abstol=abstol, # Tolerance for the absolute change of the primal-dual residuals
  maxiter=1000) # Maximum number of inner ADMM iterations

## SINGLE FIT ----

### BL fit ----
{
  time_init <- proc.time()
  fit_1run_BL <- fit_logit_splasso(y, X, D, type='BL', beta_start=beta0, 
                                   lambda=lambda, alpha=alpha, eps=eps, 
                                   intercept=intercept, maxiter=maxiter, 
                                   abstol=objtol, reltol=reltol, etatol=etatol,
                                   verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_1run_BL$exetime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  time_init <- proc.time()
  fit_1run_PG <- fit_logit_splasso(y, X, D, type='PG', beta_start=beta0, 
                                   lambda=lambda, alpha=alpha, eps=eps, 
                                   intercept=intercept, maxiter=maxiter, 
                                   abstol=objtol, reltol=reltol, etatol=etatol, 
                                   verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_1run_PG$exetime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  time_init <- proc.time()
  fit_1run_PQ <- fit_logit_splasso(y, X, D, type='PQ', beta_start=beta0, 
                                   lambda=lambda, alpha=alpha, eps=eps, 
                                   intercept=intercept, phi=phi, approx=approx, 
                                   maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                   etatol=etatol, verbose=verbose, 
                                   freq=freq, ctr_admm=ctr)
  fit_1run_PQ$exetime <- (proc.time() - time_init)[3]
}

### NR fit ----
{
  time_init <- proc.time()
  fit_1run_NR <- fit_logit_splasso(y, X, D, type='NR', beta_start=beta0, 
                                   lambda=lambda, alpha=alpha, eps=eps, 
                                   intercept=intercept, maxiter=maxiter, 
                                   abstol=objtol, reltol=reltol, etatol=etatol, 
                                   verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_1run_NR$exetime <- (proc.time() - time_init)[3]
}

### Summary ----

cat("BL:", fit_1run_BL$exetime, "\n",
    "PG:", fit_1run_PG$exetime, "\n",
    "PQ:", fit_1run_PQ$exetime, "\n",
    "NR:", fit_1run_NR$exetime, "\n")

# fit_1run_list <- list("BL"=fit_1run_BL, 
#                       "PG"=fit_1run_PG, 
#                       "PQ"=fit_1run_PQ, 
#                       "NR"=fit_1run_NR)

fit_1run_list <- list("BL"=fit_1run_BL, 
                      "PG"=fit_1run_PG, 
                      "PQ"=fit_1run_PQ)

df_1run_summary <- data.frame(
  method = names(fit_1run_list),
  niter = sapply(fit_1run_list, \(.) .$niter),
  exetime = sapply(fit_1run_list, \(.) .$exetime),
  timegain = sapply(fit_1run_list, \(.) {1-fit_1run_PQ$exetime/.$exetime}),
  loglik = sapply(fit_1run_list, \(.) .$loglik),
  row.names = seq(length(fit_1run_list)))

print(df_1run_summary)

if (SAVE) {
  filename <- paste0("portland_lasso_1run_alpha", lalpha, "_summary.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste0("portland_lasso_1run_alpha", lalpha, "_loglik.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  # Log-likelihood
  niter <- max(sapply(fit_1run_list, \(.) .$niter))
  objval <- matrix(NA, nrow=niter, ncol=length(fit_1run_list))
  colnames(objval) <- names(fit_1run_list)
  for (k in 1:length(fit_1run_list)) {
    objval[1:fit_1run_list[[k]]$niter, k] <- fit_1run_list[[k]]$trace$objval
  }
  matplot(objval[2:min(50, niter),], type="b", lty=1, col=COLORS, pch=MARKERS, xlab="", ylab="")
  title(xlab="Iteration", ylab="Log-Likelihood", main="Penalized Log-Likelihood")
  legend("bottomright", col=COLORS, pch=MARKERS, legend=colnames(objval))
  # Iterations
  with(df_1run_summary, {
    barplt <- barplot(niter, names.arg=method, col=COLORS, border=COLORS, las=2, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iterations", main="Number of Iterations")
  })
  # Exetime
  with(df_1run_summary, {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border=COLORS, las=2, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}

if (SAVE) {
  filename <- paste0("portland_lasso_1run_alpha", lalpha, "_map_pr.pdf")
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

ctr <- set_ctr_admm(rho=sqrt(p/n), gamma=0.01, smw=FALSE, precondition=FALSE, 
                    objtol=.1*objtol, reltol=reltol, abstol=abstol, maxiter=1000)

### BL fit ----
{
  time_init <- proc.time()
  fit_path_BL <- fit_logit_splasso_path(y, X, D, type='BL', 
                                        beta_start=beta0, lambda=lambdas, 
                                        alpha=alpha, gamma=gamma, maxiter=maxiter, 
                                        abstol=objtol, reltol=reltol, etatol=etatol, 
                                        verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_path_BL$tottime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  time_init <- proc.time()
  fit_path_PG <- fit_logit_splasso_path(y, X, D, type='PG', 
                                        beta_start=beta0, lambda=lambdas, 
                                        alpha=alpha, gamma=gamma, maxiter=maxiter, 
                                        abstol=objtol, reltol=reltol, etatol=etatol, 
                                        verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_path_PG$tottime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  time_init <- proc.time()
  fit_path_PQ <- fit_logit_splasso_path(y, X, D, type='PQ', 
                                        beta_start=beta0, lambda=lambdas, 
                                        alpha=alpha, gamma=gamma, maxiter=maxiter, 
                                        abstol=objtol, reltol=reltol, etatol=etatol, 
                                        verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_path_PQ$tottime <- (proc.time() - time_init)[3]
}

### NR fit ----
{
  time_init <- proc.time()
  fit_path_NR <- fit_logit_splasso_path(y, X, D, type='NR', 
                                        beta_start=beta0, lambda=lambdas, 
                                        alpha=alpha, gamma=1., maxiter=maxiter, 
                                        abstol=objtol, reltol=reltol, etatol=etatol, 
                                        verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_path_NR$tottime <- (proc.time() - time_init)[3]
}

### Summary ----

cat("BL:", fit_path_BL$tottime, "\n",
    "PG:", fit_path_PG$tottime, "\n",
    "PQ:", fit_path_PQ$tottime, "\n",
    "NR:", fit_path_NR$tottime, "\n")

# fit_path_list <- list("BL"=fit_path_BL, 
#                       "PG"=fit_path_PG, 
#                       "PQ"=fit_path_PQ, 
#                       "NR"=fit_path_NR)

fit_path_list <- list("BL"=fit_path_BL, 
                      "PG"=fit_path_PG, 
                      "PQ"=fit_path_PQ)

df_path_summary <- data.frame(
  method = names(fit_path_list),
  niter = sapply(fit_path_list, \(.) sum(.$niter)),
  exetime = sapply(fit_path_list, \(.) .$tottime),
  timegain = sapply(fit_path_list, \(.) {1-fit_path_PQ$tottime/.$tottime}),
  loglik = sapply(fit_path_list, \(.) mean(.$loglik)),
  row.names = seq(length(fit_path_list)))

print(df_path_summary)

if (SAVE) {
  filename <- paste("portland_lasso_path_alpha", lalpha, "_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_summary, file=filepath, row.names=FALSE)
}

df_path_extended <- data.frame(
  method = rep(names(fit_path_list), each=length(lambdas)),
  alpha = rep(alpha, times=length(fit_path_list)*length(lambdas)),
  lambda = c(sapply(fit_path_list, \(.) .$lambda)),
  niter = c(sapply(fit_path_list, \(.) .$niter)),
  exetime = c(sapply(fit_path_list, \(.) .$exetime)),
  timegain = c(sapply(fit_path_list, \(.) {1-fit_path_PQ$exetime/.$exetime})),
  loglik = c(sapply(fit_path_list, \(.) .$loglik)))

print(df_path_extended)

if (SAVE) {
  filename <- paste("portland_lasso_path_alpha", lalpha, "_extended.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_extended, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste0("portland_lasso_path_alpha", lalpha, "_timegain.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(df_path_summary, {
    par(mfrow=c(1,3))
    # Total number of iterations
    barplt <- barplot(niter, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), round(niter, 2))
    title(ylab="Iterations", main="Total number of iterations")
    # Total execution time
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Total execution Time")
    # Relative time gain
    reltime <- 100*(1-exetime[3]/exetime)
    barplt <- barplot(reltime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, reltime-0.1*sign(reltime)*max(reltime), paste(floor(reltime), "%"))
    title(ylab="Time Gain", main="Time Gain")
    par(mfrow=c(1,1))
  })
  dev.off()
}

### Solution path ----

# pch <- c(15:18)
# col <- c(2:4,7)

plot_fit_path(fit_path_list, field="dev", main="Deviance", position="bottomright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="logdet", main="Hessian Log-Determinant", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="reml", main="Restricted Log-Likelihood", position="bottomright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="norm", main="Parameter Norm", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="pnorm", main="Penalty Seminorm", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="niter", main="Log Number of Iterations", position="topright", log=TRUE, pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", pch=MARKERS, col=COLORS, lty=1)
plot_fit_path(fit_path_list, field="exetime", main="Log Execution Time", position="topright", log=TRUE, pch=MARKERS, col=COLORS, lty=1)


if (SAVE) {
  filename <- paste("portland_lasso_path_alpha", lalpha, "_exetime.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", log=FALSE, pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", log=FALSE, pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="niter", main="Log Number of Iterations", position="topright", log=TRUE, pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Log Execution Time", position="topright", log=TRUE, pch=MARKERS, col=COLORS, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

if (SAVE) {
  filename <- paste("portland_lasso_path_alpha", lalpha, "_loglik.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="reml", main="Hessian Log-Determinant", position="bottomright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="norm", main="Coefficient Norm", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="pnorm", main="Penalty Seminorm", position="topright", pch=MARKERS, col=COLORS, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

### Spatial maps ----

if (SAVE) {
  filename <- paste0("portland_lasso_path_alpha", lalpha, "_map_pr.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(2,2))
    fit_path_list_tmp <- list(BL=fit_path_BL, PG=fit_path_PG, PQ=fit_path_PQ, NR=fit_path_NR)
    
    k <- sapply(fit_path_list_tmp, \(.) which.max(.$reml))
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

ctr <- set_ctr_admm(rho=sqrt(p/n), gamma=0.01, smw=FALSE, precondition=FALSE, 
                    objtol=.1*objtol, reltol=reltol, abstol=abstol, maxiter=1000)

### BL fit ----
{
  time_init <- proc.time()
  fit_cv_BL <- fit_logit_splasso_cv(y, X, D, type='BL', nfold=nfold, seed=seed, 
                                    beta_start=beta0, lambda=lambdas, alpha=alpha, 
                                    gamma=gamma, maxiter=maxiter, abstol=objtol, 
                                    reltol=reltol, etatol=etatol, verbose=verbose, 
                                    freq=freq, ctr_admm=ctr)
  fit_cv_BL$tottime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  time_init <- proc.time()
  fit_cv_PG <- fit_logit_splasso_cv(y, X, D, type='PG', nfold=nfold, seed=seed, 
                                    beta_start=beta0, lambda=lambdas, alpha=alpha, 
                                    gamma=gamma, maxiter=maxiter, abstol=objtol, 
                                    reltol=reltol, etatol=etatol, verbose=verbose, 
                                    freq=freq, ctr_admm=ctr)
  fit_cv_PG$tottime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  time_init <- proc.time()
  fit_cv_PQ <- fit_logit_splasso_cv(y, X, D, type='PQ', nfold=nfold, seed=seed, 
                                    beta_start=beta0, lambda=lambdas, alpha=alpha, 
                                    gamma=gamma, maxiter=maxiter, abstol=objtol, 
                                    reltol=reltol, etatol=etatol, verbose=verbose, 
                                    freq=freq, ctr_admm=ctr)
  fit_cv_PQ$tottime <- (proc.time() - time_init)[3]
}

### NR fit ----
{
  time_init <- proc.time()
  fit_cv_NR <- fit_logit_splasso_cv(y, X, D, type='NR', nfold=nfold, seed=seed, 
                                    beta_start=beta0, lambda=lambdas, alpha=alpha, 
                                    gamma=gamma, maxiter=maxiter, abstol=objtol, 
                                    reltol=reltol, etatol=etatol, verbose=verbose, 
                                    freq=freq, ctr_admm=ctr)
  fit_cv_NR$tottime <- (proc.time() - time_init)[3]
}

### Summary ----

cat("BL:", fit_cv_BL$tottime, "\n",
    "PG:", fit_cv_PG$tottime, "\n",
    "PQ:", fit_cv_PQ$tottime, "\n",
    "NR:", fit_cv_NR$tottime, "\n")

fit_cv_list <- list("BL"=fit_cv_BL, 
                    "PG"=fit_cv_PG, 
                    "PQ"=fit_cv_PQ, 
                    "NR"=fit_cv_NR)

df_cv_extended <- rbind(cbind(method="BL", fit_cv_BL$cv), 
                        cbind(method="PG", fit_cv_PG$cv), 
                        cbind(method="PQ", fit_cv_PQ$cv), 
                        cbind(method="NR", fit_cv_NR$cv))

print(df_cv_extended)

if (SAVE) {
  filename <- paste0("portland_lasso_cv_alpha", lalpha, "_extended.csv")
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
  filename <- paste0("portland_lasso_cv_alpha", lalpha, "_summary.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_cv_summary, file=filepath, row.names=FALSE)
}


if (SAVE) {
  filename <- paste0("portland_lasso_cv_alpha", lalpha, "_timegain.pdf")
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
  filename <- paste0("portland_lasso_cv_alpha", lalpha, "_loglik.pdf")
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
  filename <- paste0("portland_lasso_cv_alpha", lalpha, "_map_pr.pdf")
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

## FINAL SUMMARY ----

if (FALSE) {
  df_all_path_summary <- data.frame()
  for (alpha in alphas) {
    lalpha <- as.character(100*alpha)
    lalpha <- ifelse(100*alpha>10, lalpha, paste0("0", lalpha))
    filename <- paste0("portland_lasso_path_alpha", lalpha, "_summary.CSV")
    filepath <- paste(CSVPATH, filename, sep="/")
    df_tmp <- read.csv2(file=filepath)
    df_all_path_summary <- rbind(df_all_path_summary, df_tmp)
  }
  
  df_all_path_aggregated <- df_all_path_summary %>% 
    filter(method != "NR") %>% group_by(method) %>% 
    summarise(niter=29*sum(niter), exetime=sum(exetime)) %>% 
    as.data.frame()
  
  df_all_path_aggregated$timegain <- with(df_all_path_aggregated, {
    100*(1-exetime[method=="PQ"]/exetime)
  })
  
  print(df_all_path_aggregated)
}


## END OF FILE ----
