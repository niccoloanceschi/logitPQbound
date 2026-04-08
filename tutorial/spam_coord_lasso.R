
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
DATAPATH <- "data/Spam"
SAVEPATH <- "tutorial/results/Spam"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

DATALAB <- "spam"

# Preprocessing settings
RESCALE <- TRUE

# Which lambda ("Custom", "Scales", or "CV")
LAMBDA <- "Scaled"
NREP <- 1L
ALPHA <- 1.0
LALPHA <- as.character(100*ALPHA)
LALPHA <- ifelse(100*ALPHA>10, LALPHA, paste0("0", LALPHA))

# Plot settings
MARKERS <- c(15:18,1,2,5,6)
COLORS <- c(2:4,7,6,5)
METHODS <- c("BL", "PG", "PQ", "NR")

options(digits=5)

## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## DATA LOAD ----

load(paste(DATAPATH, "spam_train.RData", sep="/"))

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
# lambdas <- 10^seq(-3, +2, by=0.25) # for solution path
lambdas <- 10^seq(-3,+5, by=.25) # for solution path
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
etatol <- 1. * sqrt(p) * mean(sqrt(rowSums(X^2)))
maxiter <- 7000L

# Progress output
verbose <- FALSE
freq <- 100

# Initial values
beta_start <- c(qlogis(mean(y)), rep(0, times=p-1))

## SINGLE FIT ----

### BL fit ----
{
  cat("BL...\n") 
  time_init <- proc.time()
  fit_1run_BL <- fit_logit_enet_coord(y, X, type='BL', beta_start=beta_start, 
                                      lambda=lambda, alpha=ALPHA, eps=eps, phi=phi, 
                                      intercept=intercept, maxiter=maxiter, 
                                      abstol=objtol, reltol=reltol, etatol=etatol,
                                      verbose=verbose, freq=freq)
  fit_1run_BL$exetime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  cat("PG...\n") 
  time_init <- proc.time()
  fit_1run_PG <- fit_logit_enet_coord(y, X, type='PG', beta_start=beta_start, 
                                      lambda=lambda, alpha=ALPHA, eps=eps, phi=phi, 
                                      intercept=intercept, maxiter=maxiter, 
                                      abstol=objtol, reltol=reltol, etatol=etatol, 
                                      verbose=verbose, freq=freq)
  fit_1run_PG$exetime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  cat("PQ...\n") 
  time_init <- proc.time()
  fit_1run_PQ <- fit_logit_enet_coord(y, X, type='PQ', beta_start=beta_start, 
                                      lambda=lambda, alpha=ALPHA, eps=eps, phi=phi, 
                                      intercept=intercept, maxiter=maxiter, 
                                      abstol=objtol, reltol=reltol, etatol=etatol, 
                                      verbose=verbose, freq=freq)
  fit_1run_PQ$exetime <- (proc.time() - time_init)[3]
}

### NR fit ----
{
  cat("NR...\n") 
  time_init <- proc.time()
  fit_1run_NR <- fit_logit_enet_coord(y, X, type='NR', beta_start=beta_start, 
                                      lambda=lambda, alpha=ALPHA, eps=eps, phi=phi, 
                                      intercept=intercept, maxiter=maxiter, 
                                      abstol=objtol, reltol=reltol, etatol=etatol, 
                                      verbose=verbose, freq=freq)
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
  filename <- paste(DATALAB, "_lasso_1run_alpha", LALPHA,"_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_1run_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste(DATALAB, "_lasso_1run_alpha", LALPHA,"_fit.RDS", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  saveRDS(fit_1run_list, file=filepath)
}

if (SAVE) {
  filename <- paste(DATALAB, "_lasso_1run_alpha", LALPHA,"_loglik.pdf", sep="")
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
    barplt <- barplot(niter, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iterations", main="Number of Iterations")
  })
  # Exetime
  with(df_1run_summary, {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border=COLORS, xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}


if (SAVE) {
  filename <- paste(DATALAB, "_lasso_1run_alpha", LALPHA,"_pairs.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 7; width <- 9; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  betas <- sapply(fit_1run_list, \(.) .$beta)
  pairs(betas[-1,], lower.panel=panel_cor, upper.panel=panel_linear)
  dev.off()
}

## SOLUTION PATH ----

### BL fit ----
{
  cat("BL...")
  time_init <- proc.time()
  fit_path_BL <- fit_logit_enet_path(y, X, type='BL', beta_start=NULL, 
                                     lambda=lambdas, alpha=ALPHA, eps=eps, 
                                     phi=phi, intercept=intercept, maxiter=maxiter, 
                                     abstol=objtol, reltol=reltol, etatol=etatol,
                                     verbose=verbose, freq=freq)
  fit_path_BL$tottime <- (proc.time() - time_init)[3]
}

### PG fit ----
{
  cat("PG...")
  time_init <- proc.time()
  fit_path_PG <- fit_logit_enet_path(y, X, type='PG', beta_start=NULL, 
                                     lambda=lambdas, alpha=ALPHA, eps=eps, 
                                     phi=phi, intercept=intercept, maxiter=maxiter, 
                                     abstol=objtol, reltol=reltol, etatol=etatol, 
                                     verbose=verbose, freq=freq)
  fit_path_PG$tottime <- (proc.time() - time_init)[3]
}

### PQ fit ----
{
  cat("PQ...") 
  time_init <- proc.time()
  fit_path_PQ <- fit_logit_enet_path(y, X, type='PQ', beta_start=NULL, 
                                     lambda=lambdas, alpha=ALPHA, eps=eps, 
                                     phi=phi, intercept=intercept, maxiter=maxiter, 
                                     abstol=objtol, reltol=reltol, etatol=etatol, 
                                     verbose=verbose, freq=freq)
  fit_path_PQ$tottime <- (proc.time() - time_init)[3]
}

### NR fit ----
{
  cat("NR...")
  time_init <- proc.time()
  fit_path_NR <- fit_logit_enet_path(y, X, type='NR', beta_start=beta_start, 
                                     lambda=lambdas, alpha=ALPHA, eps=eps, 
                                     phi=phi, intercept=intercept, maxiter=maxiter, 
                                     abstol=objtol, reltol=reltol, etatol=etatol, 
                                     verbose=verbose, freq=freq)
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


if (SAVE) {
  filename <- paste(DATALAB, "_lasso_path_alpha", LALPHA,"_fit.RDS", sep="")
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
  filename <- paste0(DATALAB, "_lasso_path_alpha", LALPHA, "_summary.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_summary, file=filepath, row.names=FALSE)
}


df_path_extended <- data.frame(
  method = rep(names(fit_path_list), each=length(lambdas)),
  alpha = rep(ALPHA, times=length(fit_path_list)*length(lambdas)),
  lambda = c(sapply(fit_path_list, \(.) .$lambda)),
  niter = c(sapply(fit_path_list, \(.) .$niter)),
  exetime = c(sapply(fit_path_list, \(.) .$exetime)),
  timegain = c(sapply(fit_path_list, \(.) {1-fit_path_PQ$exetime/.$exetime})),
  loglik = c(sapply(fit_path_list, \(.) .$loglik)))

print(df_path_extended)

if (SAVE) {
  filename <- paste0(DATALAB, "_lasso_path_alpha", LALPHA, "_extended.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_extended, file=filepath, row.names=FALSE)
}


if (SAVE) {
  filename <- paste0(DATALAB, "_lasso_path_alpha", LALPHA, "_timegain.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 8; zoom <- 1
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

if (SAVE) {
  filename <- paste(DATALAB, "_lasso_path_alpha", LALPHA, "_exetime.pdf", sep="")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 6; zoom <- 2
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(2,2))
  plot_fit_path(fit_path_list, field="niter", main="Number of Iterations", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="exetime", main="Execution Time", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="loglik", main="Penalized Log-Likelihood", position="topright", pch=MARKERS, col=COLORS, lty=1)
  plot_fit_path(fit_path_list, field="pnorm", main="Coefficient Norm", position="topright", pch=MARKERS, col=COLORS, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

## END OF FILE ----




