
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

# Plotting colors and markers
COLORS <- c(2:5)
MARKERS <- c(15:18)

## PORTLAND DATA ----

load(paste(DATAPATH, "PortlandData.RData", sep="/"))
load(paste(DATAPATH, "PortlandPred.RData", sep="/"))

## UTILITY FUNCTIONS ----

source("tutorial/tutorial_utils.R")

## VB-PDE FIT ----

n <- nrow(psi)
p <- ncol(psi)
y <- .5*obs+.5
X <- psi
D <- (1/sqrt(colSums(mass))) * stiff
P <- Matrix::crossprod(D)
normD <- sqrt(sum(D^2)/p) # == tr(Dt*D) / p
lambda <- 0.00004
eps <- 1e-8
intercept <- FALSE
abstol <- 1e-7
reltol <- 1e-6
maxiter <- 1000L
verbose <- FALSE
freq <- 10L
hpar <- list(A=1.0, B=1e+05, nu=3.0)
seed <- 123456
beta0 <- NULL

### ---- BL-VB ----
{
  cat("\n BL-VB")
  time_init <- proc.time()
  fit_BL <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="BL", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_BL$exetime <- (proc.time() - time_init)[3]
}

### ---- PG-VB ----
{
  cat("\n PG-VB")
  time_init <- proc.time()
  fit_PG <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="PG", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_PG$exetime <- (proc.time() - time_init)[3]
}

### ---- PQ-VB ----
{
  cat("\n PQ-VB")
  time_init <- proc.time()
  fit_PQ <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="PQ", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_PQ$exetime <- (proc.time() - time_init)[3]
}

### ---- MCMC ----
{
  cat("\n MCMC")
  time_init <- proc.time()
  fit_MC <- fit_logit_mcmc_adj(y=y, X=X, D=D, solver="sparse",
                               lambda=lambda, eps=eps, hpar=hpar, 
                               intercept=intercept, scaling=TRUE, nscaling=3,
                               maxiter=20000L, verbose=TRUE, freq=500L, seed=seed)
  fit_MC$exetime <- (proc.time() - time_init)[3]
}

## RESULTS ----

fit_vb_coeff <- array(NA, dim = c(p, 2, 4))
fit_vb_coeff[ ,1, 1] <- fit_BL$summary$beta$mean
fit_vb_coeff[ ,1, 2] <- fit_PG$summary$beta$mean
fit_vb_coeff[ ,1, 3] <- fit_PQ$summary$beta$mean
fit_vb_coeff[ ,1, 4] <- fit_MC$summary$beta$mean
fit_vb_coeff[ ,2, 1] <- diag(fit_BL$summary$beta$var)
fit_vb_coeff[ ,2, 2] <- diag(fit_PG$summary$beta$var)
fit_vb_coeff[ ,2, 3] <- diag(fit_PQ$summary$beta$var)
fit_vb_coeff[ ,2, 4] <- fit_MC$summary$beta$var
dimnames(fit_vb_coeff) <- list(beta = c(1:p),
                               par = c("mean", "var"),
                               method = c("BL-VB", "PG-VB", "PQ-VB", "MCMC"))

if (SAVE) {
  filename <- paste(DATALAB, "_mfvb_coeff.RData", sep="")
  filepath <- paste(RDSPATH, filename, sep="/")
  save(fit_vb_coeff, file=filepath)
}

### Summary measures ----

# Number of iterations
niter_BL <- fit_BL$state$niter
niter_PG <- fit_PG$state$niter 
niter_PQ <- fit_PQ$state$niter 
niter_MC <- fit_MC$call$maxiter

# Execution times
exetime_BL <- fit_BL$exetime
exetime_PG <- fit_PG$exetime
exetime_PQ <- fit_PQ$exetime
exetime_MC <- fit_MC$exetime

# ELBO
elbo_BL <- fit_BL$state$elbo
elbo_PG <- fit_PG$state$elbo 
elbo_PQ <- fit_PQ$state$elbo 
elbo_MC <- NaN

# Posterior means
mu_BL <- fit_BL$summary$beta$mean
mu_PG <- fit_PG$summary$beta$mean
mu_PQ <- fit_PQ$summary$beta$mean
mu_MC <- fit_MC$summary$beta$mean

# Posterior variances
var_BL <- diag(fit_BL$summary$beta$var)
var_PG <- diag(fit_PG$summary$beta$var)
var_PQ <- diag(fit_PQ$summary$beta$var)
var_MC <- fit_MC$summary$beta$var

# Accuracy scores (basis parameters)
acc_beta_BL <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_BL[j], var_BL[j])})
acc_beta_PG <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_PG[j], var_PG[j])})
acc_beta_PQ <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_PQ[j], var_PQ[j])})
acc_beta_MC <- rep(1.0, times=p)

# Accuracy scores (predicted field)
loc <- expand.grid(x=lon_new, y=lat_new)
out <- sp::point.in.polygon(loc[,1], loc[,2], nodes[,1], nodes[,2])

acc_field_BL <- as.vector(as.matrix(psi_new %*% acc_beta_BL))
acc_field_PG <- as.vector(as.matrix(psi_new %*% acc_beta_PG))
acc_field_PQ <- as.vector(as.matrix(psi_new %*% acc_beta_PQ))
acc_field_MC <- NA

acc_field_BL[out == 0] <- NA
acc_field_PG[out == 0] <- NA
acc_field_PQ[out == 0] <- NA

# Average accuracy/error scores
tvd_pdf_beta  <- 1-rowMeans(rbind(acc_beta_BL, acc_beta_PG, acc_beta_PQ, acc_beta_MC))
tvd_pdf_field <- 1-rowMeans(rbind(acc_field_BL, acc_field_PG, acc_field_PQ, acc_field_MC), na.rm=TRUE)

# Summary table
fit_vb_summary <- data.frame(
  method = c("BL-VB", "PG-VB", "PQ-VB", "MCMC"),
  niter = c(niter_BL, niter_PG, niter_PQ, niter_MC),
  exetime = round(c(exetime_BL, exetime_PG, exetime_PQ, exetime_MC), 2),
  timeratio = round(c(exetime_BL, exetime_PG, exetime_PQ, exetime_MC)/exetime_PQ, 4),
  elbo = round(c(elbo_BL, elbo_PG, elbo_PQ, elbo_MC), 4),
  tv = round(tvd_pdf_field, 6),
  row.names = c(1:4)
)

if (SHOW) {
  cat("\n Summary \n")
  cat(rep("-", 50), "\n", sep="", collapse="")
  print(fit_vb_summary)
  cat(rep("-", 50), "\n", sep="", collapse="")
}

if (SAVE) {
  filename <- "portland_mfvb_summary.csv"
  filepath <- paste(CSVPATH, filename, sep="/")
  write.csv2(fit_vb_summary, file=filepath, row.names=FALSE)
}

### summary plots ----

if (SAVE) {
  filename <- "portland_mfvb_elbo.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  layout(matrix(1:3, nrow = 1), widths = c(2, 1, 1))
  # ELBO
  xlim <- c(0, min(50, max(niter_BL, niter_PG, niter_PQ)))
  ylim <- range(c(fit_BL$trace$objval, fit_PG$trace$objval, fit_PQ$trace$objval))/n
  plot(0, 0, pch=19, col=3, xlim=xlim, ylim=ylim, type="o", xlab="", ylab="")
  points(fit_BL$trace$iter, fit_BL$trace$objval/n, pch=MARKERS[1], col=COLORS[1], type="o")
  points(fit_PG$trace$iter, fit_PG$trace$objval/n, pch=MARKERS[2], col=COLORS[2], type="o")
  points(fit_PQ$trace$iter, fit_PQ$trace$objval/n, pch=MARKERS[3], col=COLORS[3], type="o")
  title(xlab="Iterations", ylab="ELBO", main="Evidence lower bound")
  legend("bottomright", pch=c(15,19,17), col=2:4, legend=c("BL-VB", "PG-VB", "PQ-VB"))
  # n-iter
  with(fit_vb_summary[1:3,], {
    barplt <- barplot(niter, names.arg=method, col=COLORS, border="white", xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iteration", main="Number of Iterations")
  })
  # exe-time
  with(fit_vb_summary[1:3,], {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border="white", xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}

## PAPER PLOT ----

if (SAVE) {
  palette <- "inferno"
  filename <- paste0("portland_mfvb_field_tvd.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.25

  loc <- expand.grid(x=lon_new, y=lat_new)
  out <- sp::point.in.polygon(loc[,1], loc[,2], nodes[,1], nodes[,2])
  dat <- data.frame(x=locs[,1], y=locs[,2], z=obs)
  tvd <- 1-cbind("BL-VB"=acc_beta_BL, "PG-VB"=acc_beta_PG, "PQ-VB"=acc_beta_PQ)
  tvd <- as.matrix(psi_new %*% tvd)
  tvd[out == 0, ] <- NA

  plt <- ggplot_field(lon_new, lat_new, tvd, ngrid=ngrid, 
                      locs=dat, fun=NULL, palette=palette)
  plt <- plt + theme(axis.title=element_blank(), 
                     axis.ticks=element_blank(), 
                     axis.text=element_blank(), 
                     strip.text=element_text(size=18),
                     legend.title=element_text(size=14),
                     legend.text=element_text(size=14),
                     legend.key.height=unit(.125, units="npc"),
                     legend.key.width=unit(.0125, units="npc"))
  
  if (SHOW) print(plt)
  
  ggplot2::ggsave(filename=filename, plot=plt, path=IMGPATH, 
                  width=zoom*width, height=zoom*height)
}

## END OF FILE ----
