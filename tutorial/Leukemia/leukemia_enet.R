
## PACKAGE IMPORT ----

library(logitPQbound)
library(dplyr)
library(ggplot2)

## GLOBAL VARIABLES ----

# Saving and showing flags
SHOW <- TRUE
SAVE <- TRUE

# Global paths
DATALAB <- "leukemia"
DATAPATH <- "data/Leukemia"
SAVEPATH <- "tutorial/Leukemia"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

# Preprocessing settings
RESCALE <- TRUE

# E-net mixing weight
ALPHA <- 0.8
LALPHA <- as.character(100*ALPHA)
LALPHA <- ifelse(100*ALPHA>10, LALPHA, paste0("0", LALPHA))

options(digits=5)

## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## DATA LOAD ----

load(paste(DATAPATH, "Leukemia.RData", sep="/"))
str(Leukemia)

## PREPROCESSING ----

n <- nrow(Leukemia$x)
p <- ncol(Leukemia$x)+1

X <- as.matrix(Leukemia$x)
y <- as.vector(Leukemia$y)

X <- cbind(intercept=1, X)
D <- diag(p)[-1,]

if (RESCALE) {
  X[,-1] <- scale(X[,-1])
}

rm(Leukemia)
gc()

## MODEL SET-UP ----

# Penalty parameters
lambdas <- 10^seq(-3, +2, by=0.25)

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

### Summary ----

fit_path_list <- list("BL"=fit_path_BL, "PG"=fit_path_PG, "PQ"=fit_path_PQ)
fit_path_summary <- get_path_summary(fit_path_list, ALPHA)
fit_path_coeff <- get_path_coeff(fit_path_list, lambdas)

if (SHOW) {
  cat("\n Summary \n")
  cat(rep("-", 50), "\n", sep="", collapse="")
  print(fit_path_summary)
  cat(rep("-", 50), "\n", sep="", collapse="")
}

if (SAVE) {
  filename <- paste(DATALAB, "_enet", LALPHA, "_path_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(fit_path_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste(DATALAB, "_enet", LALPHA, "_path_coeff.RData", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  save(ALPHA, lambdas, fit_path_coeff, file=filepath)
}

## END OF FILE ----
