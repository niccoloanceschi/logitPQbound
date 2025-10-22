
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
SAVE <- FALSE

# Global paths
DATAPATH <- "data/Alzheimer"
SAVEPATH <- "tutorial/results/Alzheimer"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

# Preprocessing settings
RESCALE <- TRUE

# Which lambda ("Custom", "Scales", or "CV")
LAMBDA <- "Scaled"
NREP <- 10L

# Plot settings
MARKERS <- c(15:18)
COLORS <- c(2:4,7,6)
METHODS <- c("BL", "PG", "PQ", "NR")

options(digits=5)

## DATA LOAD ----

load(paste(DATAPATH, "Alzheimer.RData", sep="/"))

## PREPROCESSING ----

X <- as.matrix(X)
y <- as.vector(y)

n <- length(y)
p <- ncol(X)

idx_train <- sample(1:n, size=floor(.8*n), replace=FALSE)
idx_test <- setdiff(1:n, idx_train)

n_train <- length(idx_train)
n_test  <- length(idx_test)

X_train <- X[idx_train,]
y_train <- y[idx_train]

X_test <- X[idx_test,]
y_test <- y[idx_test]

if (RESCALE) {
  X_means <- apply(X_train[,-1], 2, mean)
  X_stds <- apply(X_train[,-1], 2, sd)
  
  X_train[,-1] <- t((t(X_train[,-1]) - X_means) / X_stds)
  X_test[,-1] <- t((t(X_test[,-1]) - X_means) / X_stds)
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
beta_start <- rep(0, times=p)

## MODEL FIT ----

df_fit_extended <- data.frame()
for(r in 1:NREP){
  cat("- replication: ", r, "/", NREP, " | ", sep="")
  
  ### BL fit
  cat(" BL...")
  exetime_BL <- proc.time()
  fit_BL <- fit_logit_ridge(y, X, type='BL', beta_start=beta_start, 
                            lambda=lambda, eps=eps, intercept=intercept, 
                            phi=phi, maxiter=maxiter, abstol=objtol, 
                            reltol=reltol, etatol=etatol,
                            verbose=verbose, 
                            freq=freq)
  exetime_BL <- (proc.time() - exetime_BL)[3]
  
  ### PG fit
  cat(" PG...")
  exetime_PG <- proc.time()
  fit_PG <- fit_logit_ridge(y, X, type='PG', beta_start=beta_start, 
                            lambda=lambda, eps=eps, intercept=intercept, 
                            phi=phi, maxiter=maxiter, abstol=objtol, 
                            reltol=reltol, etatol=etatol, 
                            verbose=verbose, freq=freq)
  exetime_PG <- (proc.time() - exetime_PG)[3]
  
  ### PQ fit
  cat(" PQ...")
  exetime_PQ <- proc.time()
  fit_PQ <- fit_logit_ridge(y, X, type='PQ', beta_start=beta_start, 
                            lambda=lambda, eps=eps, intercept=intercept, 
                            phi=phi, maxiter=maxiter, abstol=objtol, 
                            reltol=reltol, etatol=etatol, 
                            verbose=verbose, freq=freq)
  exetime_PQ <- (proc.time() - exetime_PQ)[3]
  
  ### NR fit
  cat(" NR...\n")
  exetime_NR <- proc.time()
  fit_NR <- fit_logit_ridge(y, X, type='NR', beta_start=beta_start, 
                            lambda=lambda, eps=eps, intercept=intercept, 
                            phi=phi, maxiter=maxiter, abstol=objtol, 
                            reltol=reltol, etatol=etatol, 
                            verbose=verbose, freq=freq)
  exetime_NR <- (proc.time() - exetime_NR)[3]
  
  ### Summary
  df_fit_tmp <- data.frame(
    replication = rep(r, times=4),
    method = c("BL", "PG", "PQ", "NR"),
    niter = c(fit_BL$niter, fit_PG$niter, fit_PQ$niter, fit_NR$niter),
    exetime = c(exetime_BL, exetime_PG, exetime_PQ, exetime_NR),
    loglik = c(fit_BL$loglik, fit_PG$loglik, fit_PQ$loglik, fit_NR$loglik))
  
  df_fit_extended <- rbind(df_fit_extended, df_fit_tmp)
  rm(df_fit_tmp); gc()
}

## SUMMARY ----

df_fit_summary <- 
  df_fit_extended %>% 
  dplyr::mutate(method=factor(method, levels=c("BL", "PG", "PQ", "NR"))) %>% 
  dplyr::group_by(method) %>% 
  dplyr::summarise(niter=mean(niter), 
                   exetime=mean(exetime), 
                   loglik=mean(loglik)) %>% 
  as.data.frame()

print(df_fit_summary)

if (SAVE) {
  filename <- "alzheimer_ridge_extended.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_fit_extended, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- "alzheimer_ridge_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_fit_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  fit_list <- list(BL=fit_BL, PG=fit_PG, PQ=fit_PQ, NR=fit_NR)
  filename <- "alzheimer_ridge_fit.RDS"
  filepath <- paste(RDSPATH, filename, sep="/")
  saveRDS(fit_list, file=filepath)
}

if (SAVE) {
  filename <- "alzheimer_ridge_loglik.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  # Log-likelihood
  niter <- max(sapply(fit_list, \(.) .$niter))
  objval <- matrix(NA, nrow=niter, ncol=4)
  objval[1:fit_BL$niter, 1] <- fit_BL$trace$objval
  objval[1:fit_PG$niter, 2] <- fit_PG$trace$objval
  objval[1:fit_PQ$niter, 3] <- fit_PQ$trace$objval
  objval[1:fit_NR$niter, 4] <- fit_NR$trace$objval
  colnames(objval) <- METHODS
  matplot(objval[2:min(50, niter),], type="b", lty=1, col=COLORS, pch=MARKERS, xlab="", ylab="")
  title(xlab="Iteration", ylab="Log-Likelihood", main="Penalized Log-Likelihood")
  legend("bottomright", col=COLORS, pch=MARKERS, legend=METHODS)
  # Iterations
  with(df_fit_summary, {
    barplt <- barplot(niter, names.arg=method, col=COLORS, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iterations", main="Number of Iterations")
  })
  # Exetime
  with(df_fit_summary, {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}

