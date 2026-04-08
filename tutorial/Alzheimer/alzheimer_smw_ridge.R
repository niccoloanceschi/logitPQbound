
## PACKAGE IMPORT ----

suppressWarnings({
  # PQ-bound
  devtools::load_all()
  # library(logitPQbound)
  
  # Coordinate ascent for GLMs 
  library(glmnet)
  
  # Tidyverse
  library(dplyr)
  library(ggplot2)
  
  # Viridis colours
  library(viridis)
  
  # Sparse linear algebra
  library(Matrix)
  library(sparseinv)
  
  # Benchmarking
  library(microbenchmark)
})


## GLOBAL VARIABLES ----

# Saving and showing flags
SHOW <- TRUE
SAVE <- TRUE

# Global paths
DATAPATH <- "data/Alzheimer"
SAVEPATH <- "tutorial/Alzheimer"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

DATALAB <- "alzheimer"

# Preprocessing settings
RESCALE <- TRUE

# Which lambda ("Custom", "Scales", or "CV")
LAMBDA <- "Scaled"
NREP <- 10L
ALPHA <- 0

# Plot settings
MARKERS <- c(15:19)
COLORS <- c(2:4,7,6)

options(digits=5)


## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## DATA LOAD ----

load(paste(DATAPATH, "Alzheimer.RData", sep="/"))

## PREPROCESSING ----

n <- nrow(X)
p <- ncol(X)

y <- as.vector(y)
X <- as.matrix(X)
D <- diag(p)

if (RESCALE) {
  X[,-1] <- scale(X[,-1])
}

## MODEL SET-UP ----

# Penalty parameters
lambdas <- 10^seq(-3, +2, by=0.25) # for solution path
lambda <- .1 # for single fit

# Cross-validation Seed and n of folds
seed <- 123456
nfold <- 5

# Intercept penalty
eps <- 1e-8
intercept <- TRUE

# EDF inflation factor (for GCV computation only)
gamma <- 1.

# PQ update
phi <- 0.9

# Convergence tolerance 
objtol <- 1e-7
reltol <- 1e-2
abstol <- 1e-3
etatol <- 1.
maxiter <- 5000L

# Progress output
verbose <- FALSE
freq <- 100

# Initial values
beta_start <- c(qlogis(mean(y)), rep(0, times=p-1))

## SINGLE FIT ----

### BL fit ----
{
  cat(" BL...")
  time_init <- proc.time()
  fit_1run_BL <- fit_logit_ridge(y, X, type='BL', beta_start=beta_start, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol,
                                 verbose=verbose, freq=freq)
  fit_1run_BL$exetime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  cat(" PG...")
  time_init <- proc.time()
  fit_1run_PG <- fit_logit_ridge(y, X, type='PG', beta_start=beta_start, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol, 
                                 verbose=verbose, freq=freq)
  fit_1run_PG$exetime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  cat(" PQ...")
  time_init <- proc.time()
  fit_1run_PQ <- fit_logit_ridge(y, X, type='PQ', beta_start=beta_start, 
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

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_1run_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_1run_fit.RDS", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  saveRDS(fit_1run_list, file=filepath)
}

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_1run_loglik.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_1run_loglik(fit_1run_list, df_1run_summary, COLORS, MARKERS)
  dev.off()
}

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_1run_pairs.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 7; width <- 9; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_1run_pairs(fit_1run_list)
  dev.off()
}


## SOLUTION PATH ----

### BL fit ----
{
  cat(" BL...")
  time_init <- proc.time()
  fit_path_BL <- fit_logit_ridge_path(y, X, type='BL', beta_start=beta_start, 
                                      lambda=lambdas, eps=eps, intercept=intercept, 
                                      phi=phi, maxiter=maxiter, abstol=objtol, 
                                      reltol=reltol, etatol=etatol,
                                      verbose=verbose, freq=freq)
  fit_path_BL$tottime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  cat(" PG...")
  time_init <- proc.time()
  fit_path_PG <- fit_logit_ridge_path(y, X, type='PG', beta_start=beta_start, 
                                      lambda=lambdas, eps=eps, intercept=intercept, 
                                      phi=phi, maxiter=maxiter, abstol=objtol, 
                                      reltol=reltol, etatol=etatol, 
                                      verbose=verbose, freq=freq)
  fit_path_PG$tottime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  cat(" PQ...")
  time_init <- proc.time()
  fit_path_PQ <- fit_logit_ridge_path(y, X, type='PQ', beta_start=beta_start, 
                                      lambda=lambdas, eps=eps, intercept=intercept, 
                                      phi=phi, maxiter=maxiter, abstol=objtol, 
                                      reltol=reltol, etatol=etatol, 
                                      verbose=verbose, freq=freq)
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
  alpha = rep(ALPHA, length(fit_path_list)),
  method = names(fit_path_list),
  niter = sapply(fit_path_list, \(.) sum(.$niter)),
  exetime = sapply(fit_path_list, \(.) .$tottime),
  timegain = sapply(fit_path_list, \(.) {1-fit_path_PQ$tottime/.$tottime}),
  loglik = sapply(fit_path_list, \(.) mean(.$loglik)),
  row.names = seq(length(fit_path_list)))

print(df_path_summary)


if (SAVE) {
  filename <- paste(DATALAB, "_ridge_path_summary.csv", sep="")
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

if (SAVE) {
  filename <- paste0(DATALAB, "_ridge_path_extended.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_extended, file=filepath, row.names=FALSE)
}


if (SAVE) {
  filename <- paste0(DATALAB, "_ridge_path_timegain.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 8; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_path_timegain(df_path_summary, COLORS, MARKERS)
  dev.off()
}


### Solution path ----

if (SAVE) {
  filename <- paste(DATALAB, "_ridge_path_exetime.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="norm", main="Coefficient Norm", position="topright", pch=MARKERS, col=COLORS, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

## END OF FILE ----
