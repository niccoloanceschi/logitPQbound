
## PACKAGE IMPORT ----

suppressWarnings({
  # PQ-bound
  library(logitPQbound)
  
  # Tidyverse
  library(dplyr)
  library(ggplot2)
})


## GLOBAL VARIABLES ----

# Saving and showing flags
SHOW <- TRUE
SAVE <- TRUE

# Global paths
DATAPATH <- "data/Spam"
SAVEPATH <- "tutorial/Spam"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

DATALAB <- "spam"

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
  filename <- paste(DATALAB, "_logit_1run_fit.RDS", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  saveRDS(fit_1run_list, file=filepath)
}

if (SAVE) {
  filename <- paste(DATALAB, "_logit_1run_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}

if (FALSE) {
  filename <- paste(DATALAB, "_logit_1run_loglik.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_1run_loglik(fit_1run_list, df_1run_summary, COLORS, MARKERS)
  dev.off()
}

if (FALSE) {
  filename <- paste(DATALAB, "_logit_1run_pairs.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 7; width <- 9; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  plot_1run_pairs(fit_1run_list)
  dev.off()
}

## END OF FILE ----
