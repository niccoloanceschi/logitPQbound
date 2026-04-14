
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

load(paste(DATAPATH, "Alzheimer.RData", sep="/"))

## PREPROCESSING ----

n <- nrow(X)
p <- ncol(X)

y <- as.vector(y)
X <- as.matrix(X)
D <- diag(p)[-1,]

if (RESCALE) {
  X[,-1] <- scale(X[,-1])
}

## MODEL SET-UP ----

# Penalty parameters
lambdas <- 10^seq(-3,+2, by=.25) # for solution path
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
maxiter <- 5000L

# Progress output
verbose <- FALSE
freq <- 100

# Initial values
beta_start <- c(qlogis(mean(y)), rep(0, times=p-1))

## SOLUTION PATH ----

### BL fit ----
{
  cat("BL...\n")
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
  cat("PG...\n")
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
  cat("PQ...\n") 
  time_init <- proc.time()
  fit_path_PQ <- fit_logit_enet_path(y, X, type='PQ', beta_start=NULL, 
                                     lambda=lambdas, alpha=ALPHA, eps=eps, 
                                     phi=phi, intercept=intercept, maxiter=maxiter, 
                                     abstol=objtol, reltol=reltol, etatol=etatol, 
                                     verbose=verbose, freq=freq)
  fit_path_PQ$tottime <- (proc.time() - time_init)[3]
}

### Summary ----

fit_path_list <- list("BL"=fit_path_BL, 
                      "PG"=fit_path_PG, 
                      "PQ"=fit_path_PQ)

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
  filename <- paste0(DATALAB, "_lasso_path_summary.csv")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(df_path_summary, file=filepath, row.names=FALSE)
}

fit_path_coeff <- array(NA, dim = c(p, length(lambdas), 3))
fit_path_coeff[,,1] <- fit_path_BL$beta
fit_path_coeff[,,2] <- fit_path_PG$beta
fit_path_coeff[,,3] <- fit_path_PQ$beta
dimnames(fit_path_coeff) <- list(beta = 1:p,
                                 lambda = 1:length(lambdas),
                                 method = c("BL", "PG", "PQ"))

if (SAVE) {
  filename <- paste(DATALAB, "_lasso_path_coeff.RData", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  save(alpha, lambdas, fit_path_coeff, file=filepath)
}


## END OF FILE ----




