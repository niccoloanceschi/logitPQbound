
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
alphas <- seq(.05, 1, by=0.15) # for solution path

lambda <- .0001 # for single fit
alpha <- alphas[7] # Mixing weight of L1 and L2 penalties

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
  rho=sqrt(p/n), # Augmented Lagrangian penalty parameter
  gamma=0.01, # Proximal penalty parameter
  smw=FALSE, # Sherman-Morrison-Woodbury update
  precondition=FALSE, # Preconditioned ADMM
  objtol=.1*objtol, # Tolerance for the objective function
  reltol=reltol, # Tolerance for the relative change of the primal-dual residuals
  abstol=abstol, # Tolerance for the absolute change of the primal-dual residuals
  maxiter=1000) # Maximum number of inner ADMM iterations

## SOLUTION PATH ----

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

### Summary ----

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
  filename <- paste("portland_lasso_path_summary.csv", sep="")
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
