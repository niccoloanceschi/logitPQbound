
## PACKAGE IMPORT ----

library(logitPQbound)
library(dplyr)
library(ggplot2)
library(Matrix)
library(sparseinv)

## GLOBAL VARIABLES ----

# Saving and showing flags
SHOW <- TRUE
SAVE <- TRUE

# Global paths
DATALAB <- "portland"
DATAPATH <- "data/Portland"
SAVEPATH <- "tutorial/Portland"
RDSPATH <- paste(SAVEPATH, "rds", sep="/")
CSVPATH <- paste(SAVEPATH, "csv", sep="/")
IMGPATH <- paste(SAVEPATH, "img", sep="/")

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

# Random seed
seed <- 123456

# Intercept penalty
eps <- 1e-8
intercept <- FALSE

# PQ proximal update
phi <- 0.5

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

fit_path_BL <- list()
fit_path_PG <- list()
fit_path_PQ <- list()

df_path_extended <- data.frame()

for (k in 1:length(alphas)) {
  
  ### Set the mixing weight
  alpha <- alphas[k]
  
  ### BL fit ----
  cat("\n BL | alpha =", alpha, "\n")
  cat(rep("-", 65), "\n", sep="", collapse="")
  time_init <- proc.time()
  fit_tmp_BL <- fit_logit_splasso_path(y, X, D, type='BL', 
                                       beta_start=beta0, lambda=lambdas, 
                                       alpha=alpha, maxiter=maxiter, 
                                       abstol=objtol, reltol=reltol, etatol=etatol, 
                                       verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_tmp_BL$tottime <- (proc.time() - time_init)[3]
  fit_path_BL[[k]] <- fit_tmp_BL[-c(7:14)]
  
  ## PG fit ----
  cat("\n PG | alpha =", alpha, "\n")
  cat(rep("-", 65), "\n", sep="", collapse="")
  time_init <- proc.time()
  fit_tmp_PG <- fit_logit_splasso_path(y, X, D, type='PG', 
                                       beta_start=beta0, lambda=lambdas, 
                                       alpha=alpha, maxiter=maxiter, 
                                       abstol=objtol, reltol=reltol, etatol=etatol, 
                                       verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_tmp_PG$tottime <- (proc.time() - time_init)[3]
  fit_path_PG[[k]] <- fit_tmp_PG[-c(7:14)]
  
  ## PQ fit ----
  cat("\n PQ | alpha =", alpha, "\n")
  cat(rep("-", 65), "\n", sep="", collapse="")
  time_init <- proc.time()
  fit_tmp_PQ <- fit_logit_splasso_path(y, X, D, type='PQ', 
                                       beta_start=beta0, lambda=lambdas, 
                                       alpha=alpha, maxiter=maxiter, 
                                       abstol=objtol, reltol=reltol, etatol=etatol, 
                                       verbose=verbose, freq=freq, ctr_admm=ctr)
  fit_tmp_PQ$tottime <- (proc.time() - time_init)[3]
  fit_path_PQ[[k]] <- fit_tmp_PQ[-c(7:14)]
  
  ## Storing ----
  fit_path_tmp <- list(
    "BL"=fit_tmp_BL[-c(7:14)], 
    "PG"=fit_tmp_PG[-c(7:14)], 
    "PQ"=fit_tmp_PQ[-c(7:14)])
  
  df_path_tmp <- data.frame(
    method = names(fit_path_tmp),
    alpha = sapply(fit_path_tmp, \(.) .$alpha),
    niter = sapply(fit_path_tmp, \(.) sum(.$niter)),
    exetime = sapply(fit_path_tmp, \(.) .$tottime),
    timeratio = sapply(fit_path_tmp, \(.) {.$tottimefit_tmp_PQ$tottime}))
    # loglik = sapply(fit_path_tmp, \(.) mean(.$loglik)))
  
  ### Store the results
  df_path_extended <- rbind(df_path_extended, df_path_tmp)
  
  ### Free the memory 
  rm(fit_tmp_BL, fit_tmp_PG, fit_tmp_PQ, fit_path_tmp, df_path_tmp)
  gc()
}

### Summary ----

# Summary statistics
fit_path_summary <- df_path_extended %>% 
  group_by(method) %>% 
  summarise(niter = round(sum(niter), 4), 
            exetime = round(sum(exetime), 4),
            .groups = "drop") %>%
  mutate(time_PQ = exetime[method == "PQ"],
         timeratio = round(exetime/time_PQ, 4)) %>% 
  select(-time_PQ) %>% 
  as.data.frame()

# Estimated coefficients' array
fit_path_coeff <- array(NA, dim = c(p, length(lambdas), length(alphas), 3))
for (k in 1:length(alphas)) {
  fit_path_coeff[,,k,1] <- fit_path_BL[[k]]$beta
  fit_path_coeff[,,k,2] <- fit_path_PG[[k]]$beta
  fit_path_coeff[,,k,3] <- fit_path_PQ[[k]]$beta
}
dimnames(fit_path_coeff) <- list(
  beta = 1:p,
  lambda = 1:length(lambdas),
  alpha = 1:length(alphas),
  method = c("BL", "PG", "PQ"))

if (SHOW) {
  cat("\n Summary \n")
  cat(rep("-", 50), "\n", sep="", collapse="")
  print(fit_path_summary)
  cat(rep("-", 50), "\n", sep="", collapse="")
}

if (SAVE) {
  filename <- paste(DATALAB, "_enet2d_path_summary.csv", sep="")
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(fit_path_summary, file=filepath, row.names=FALSE)
}

if (SAVE) {
  filename <- paste(DATALAB, "_enet2d_path_coeff.RData", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  save(alphas, lambdas, fit_path_coeff, file=filepath)
}


## END OF FILE ----
