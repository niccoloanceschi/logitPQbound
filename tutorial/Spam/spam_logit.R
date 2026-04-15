
## PACKAGE IMPORT ----

library(logitPQbound)
library(dplyr)
library(ggplot2)

## GLOBAL VARIABLES ----

# Saving and showing flags
SHOW <- TRUE
SAVE <- TRUE

# Global paths
DATALAB <- "spam"
DATAPATH <- "data/Spam"
SAVEPATH <- "tutorial/Spam"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

# Preprocessing settings
RESCALE <- TRUE

# E-net mixing weight
ALPHA <- 0

options(digits=5)


## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## DATA LOAD ----

load(paste(DATAPATH, "Spam.RData", sep="/"))

## PREPROCESSING ----

n <- nrow(Spam$x)
p <- ncol(Spam$x)+1

X <- as.matrix(Spam$x)
y <- as.vector(Spam$y)

X <- cbind(intercept=1, X)
D <- diag(p)[-1,]

if (RESCALE) {
  for (j in 2:p) {
    X[,j] <- log1p(X[,j])
    mex_xj <- median(X[,j])
    iqr_xj <- IQR(X[,j])
    X[,j] <- (X[,j] - mex_xj)
    if (iqr_xj>0) X[,j] / iqr_xj
  }
}

rm(Spam)
gc()


## MODEL SET-UP ----

# Penalty parameters
lambda <- 1e-16

# Random seed (for random scan coordinate ascent)
seed <- 123456

# Intercept penalty
eps <- 1e-8
intercept <- TRUE

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
  cat("\n BL bound")
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
  cat("\n PG bound")
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
  cat("\n PQ bound\n")
  time_init <- proc.time()
  fit_1run_PQ <- fit_logit_ridge(y, X, type='PQ', beta_start=beta_start, 
                                 lambda=lambda, eps=eps, intercept=intercept, 
                                 phi=phi, maxiter=maxiter, abstol=objtol, 
                                 reltol=reltol, etatol=etatol, 
                                 verbose=verbose, freq=freq)
  fit_1run_PQ$exetime <- (proc.time() - time_init)[3]
}

### Summary ----

fit_1run_list <- list("BL"=fit_1run_BL, "PG"=fit_1run_PG, "PQ"=fit_1run_PQ)
fit_1run_summary <- get_1run_summary(fit_1run_list, ALPHA)
fit_1run_coeff <- get_1run_coeff(fit_1run_list, lambdas)


if (SHOW) {
  cat("\n Summary \n")
  cat(rep("-", 50), "\n", sep="", collapse="")
  print(fit_1run_summary)
  cat(rep("-", 50), "\n", sep="", collapse="")
}

if (SAVE) {
  filename <- paste(DATALAB, "_logit_1run_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(fit_1run_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste(DATALAB, "_logit_1run_coeff.RData", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  save(fit_1run_coeff, file=filepath)
}

## END OF FILE ----
