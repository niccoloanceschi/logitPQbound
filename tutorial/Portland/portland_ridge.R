
## PACKAGE IMPORT ----

library(logitPQbound)
library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

SHOW <- TRUE
SAVE <- TRUE

DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/Portland"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

DATALAB <- "portland"

COLORS  <- c(2:4,7)
MARKERS <- c(15:18)

## PORTLAND DATA ----

load(paste(DATAPATH, "PortlandData.RData", sep="/"))

## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## OPTIM SET-UP ----

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

# Set the random seed
set.seed(seed)

# Initial values
beta0 <- rnorm(p, mean=0, sd=sqrt(1/p))

## SINGLE FIT ----

### BL fit ----
{
  time_init <- proc.time()
  fit_1run_BL <- fit_logit_spridge(y, X, D, type='BL', beta_start=beta0, 
                                   lambda=lambda, eps=eps, intercept=intercept, 
                                   phi=phi, maxiter=maxiter, abstol=objtol, 
                                   reltol=reltol, etatol=etatol,
                                   verbose=verbose, freq=freq)
  fit_1run_BL$exetime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  time_init <- proc.time()
  fit_1run_PG <- fit_logit_spridge(y, X, D, type='PG', beta_start=beta0, 
                                   lambda=lambda, eps=eps, intercept=intercept, 
                                   phi=phi, maxiter=maxiter, abstol=objtol, 
                                   reltol=reltol, etatol=etatol, 
                                   verbose=verbose, freq=freq)
  fit_1run_PG$exetime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  time_init <- proc.time()
  fit_1run_PQ <- fit_logit_spridge(y, X, D, type='PQ', beta_start=beta0, 
                                   lambda=lambda, eps=eps, intercept=intercept, 
                                   phi=phi, maxiter=maxiter, abstol=objtol, 
                                   reltol=reltol, etatol=etatol, 
                                   verbose=verbose, freq=freq)
  fit_1run_PQ$exetime <- (proc.time() - time_init)[3]
}

### Summary ----

cat("BL:", fit_1run_BL$exetime, "\n",
    "PG:", fit_1run_PG$exetime, "\n",
    "PQ:", fit_1run_PQ$exetime, "\n")

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

if (FALSE) {
  filename <- "portland_ridge_1run_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}


if (FALSE) {
  filename <- "portland_ridge_1run_loglik.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_1run_loglik(fit_1run_list, df_1run_summary, COLORS, MARKERS)
  dev.off()
}

if (FALSE) {
  filename <- "portland_ridge_1run_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4.5; width <- 15; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(1,3))
    pred_1run_BL <- as.vector(1 / (1 + exp(- X %*% fit_1run_BL$beta)))
    pred_1run_PG <- as.vector(1 / (1 + exp(- X %*% fit_1run_PG$beta)))
    pred_1run_PQ <- as.vector(1 / (1 + exp(- X %*% fit_1run_PQ$beta)))
    
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_BL, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="BL"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_PG, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PG"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_1run_PQ, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PQ"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

## SOLUTION PATH ----

### BL fit ----
{
  time_init <- proc.time()
  fit_path_BL <- fit_logit_spridge_path(y, X, D, type='BL', beta_start=beta0, 
                                        lambda=lambdas, gamma=gamma, phi=phi, 
                                        maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                        etatol=etatol, verbose=verbose, freq=freq)
  fit_path_BL$tottime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  time_init <- proc.time()
  fit_path_PG <- fit_logit_spridge_path(y, X, D, type='PG', beta_start=beta0, 
                                        lambda=lambdas, gamma=gamma, phi=phi, 
                                        maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                        etatol=etatol, verbose=verbose, freq=freq)
  fit_path_PG$tottime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  time_init <- proc.time()
  fit_path_PQ <- fit_logit_spridge_path(y, X, D, type='PQ', beta_start=beta0, 
                                        lambda=lambdas, gamma=gamma, phi=phi, 
                                        maxiter=maxiter, abstol=objtol, reltol=reltol, 
                                        etatol=etatol, verbose=verbose, freq=freq)
  fit_path_PQ$tottime <- (proc.time() - time_init)[3]
}

### Summary ----

cat("BL:", fit_path_BL$tottime, "\n",
    "PG:", fit_path_PG$tottime, "\n",
    "PQ:", fit_path_PQ$tottime, "\n")

fit_path_list <- list("BL"=fit_path_BL, 
                      "PG"=fit_path_PG, 
                      "PQ"=fit_path_PQ)

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_path_fit.RDS", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  saveRDS(fit_path_list, file=filepath)
}

df_path_summary <- data.frame(
  method = names(fit_path_list),
  niter = sapply(fit_path_list, \(.) sum(.$niter)),
  exetime = sapply(fit_path_list, \(.) .$tottime),
  timegain = sapply(fit_path_list, \(.) {1-fit_path_PQ$tottime/.$tottime}),
  loglik = sapply(fit_path_list, \(.) mean(.$loglik)),
  row.names = seq(length(fit_path_list)))

print(df_path_summary)

if (SAVE) {
  filename <- "portland_ridge_path_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_summary, file=filepath, row.names=FALSE)
}

df_path_extended <- data.frame(
  method = rep(names(fit_path_list), each=length(lambdas)),
  alpha = rep(.0, times=length(fit_path_list)*length(lambdas)),
  lambda = c(sapply(fit_path_list, \(.) .$lambda)),
  niter = c(sapply(fit_path_list, \(.) .$niter)),
  exetime = c(sapply(fit_path_list, \(.) .$exetime)),
  timegain = c(sapply(fit_path_list, \(.) {1-fit_path_PQ$exetime/.$exetime})),
  loglik = c(sapply(fit_path_list, \(.) .$loglik)))

print(df_path_extended)

if (FALSE) {
  filename <- "portland_ridge_path_extended.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_extended, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- "portland_ridge_path_timegain.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 8; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_path_timegain(df_path_summary, COLORS, MARKERS)
  dev.off()
}

### Solution path ----

if (FALSE) {
  filename <- "portland_ridge_path_exetime.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="pnorm", main="Penalty Seminorm", position="topright", pch=MARKERS, col=COLORS, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

### Spatial maps ----

if (FALSE) {
  filename <- "portland_ridge_path_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4.5; width <- 15; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow=c(1,3))
    fit_path_list_tmp <- list(BL=fit_path_BL, PG=fit_path_PG, PQ=fit_path_PQ)
    k <- sapply(fit_path_list_tmp, \(.) which.max(.$reml))
    
    pred_path_BL <- c(plogis(as.vector(X %*% fit_path_BL$beta[,k["BL"]])))
    pred_path_PG <- c(plogis(as.vector(X %*% fit_path_PG$beta[,k["PG"]])))
    pred_path_PQ <- c(plogis(as.vector(X %*% fit_path_PQ$beta[,k["PQ"]])))
    
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_BL, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="BL"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_PG, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PG"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pred_path_PQ, pch=19, xlim=xlim, ylim=ylim, clim=clim, xlab="", ylab="", main="PQ"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

## END OF FILE ----