
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
verbose <- TRUE
freq <- 10L
hpar <- list(A=1.0, B=1e+05, nu=3.0)

# beta0 <- rnorm(p, mean=0, sd=sqrt(1/p))
beta0 <- NULL

### ---- BL-VB ----
{
  time_init <- proc.time()
  fit_BL <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="BL", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_BL$exetime <- (proc.time() - time_init)[3]
}

### ---- PG-VB ----
{
  time_init <- proc.time()
  fit_PG <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="PG", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_PG$exetime <- (proc.time() - time_init)[3]
}

### ---- PQ-VB ----
{
  time_init <- proc.time()
  fit_PQ <- fit_logit_mfvb_adj(y=y, X=X, D=D, type="PQ", beta_start=beta0, 
                               lambda=lambda, eps=eps, hpar=hpar, intercept=intercept, 
                               solver="sparse", maxiter=maxiter, abstol=abstol, 
                               reltol=reltol, verbose=verbose, freq=freq)
  fit_PQ$exetime <- (proc.time() - time_init)[3]
}

### ---- MCMC ----
{
  time_init <- proc.time()
  fit_MC <- fit_logit_mcmc_adj(y=y, X=X, D=D, solver="sparse",
                               lambda=lambda, eps=eps, hpar=hpar, 
                               intercept=intercept, scaling=TRUE, nscaling=3,
                               maxiter=20000L, verbose=verbose, freq=200L)
  fit_MC$exetime <- (proc.time() - time_init)[3]
}

## RESULTS ----

if (SAVE) {
  fit_MC_tr <- fit_MC$trace$beta
  fit_MC$trace$beta <- sqrt(rowMeans(fit_MC_tr^2))
  
  filename <- "portland_mfvb_fit.RData"
  filepath <- paste(RDSPATH, filename, sep="/")
  save(fit_BL, fit_PG, fit_PQ, fit_MC, file=filepath)
  
  fit_MC$trace$beta <- fit_MC_tr
  rm(fit_MC_tr)
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

# Posterior linear predictor means
eta_BL <- fit_BL$summary$eta$mean
eta_PG <- fit_PG$summary$eta$mean
eta_PQ <- fit_PQ$summary$eta$mean
eta_MC <- fit_MC$summary$eta$mean

trace_eta_MC <- as.matrix(tcrossprod(fit_MC$trace$beta, X))

# Posterior linear predictor variances
var_eta_BL <- fit_BL$summary$eta$var
var_eta_PG <- fit_PG$summary$eta$var
var_eta_PQ <- fit_PQ$summary$eta$var
var_eta_MC <- fit_MC$summary$eta$var

# Posterior linear predictor means
pr_BL <- 1 / (1 + exp(-eta_BL))
pr_PG <- 1 / (1 + exp(-eta_PG))
pr_PQ <- 1 / (1 + exp(-eta_PQ))
pr_MC <- 1 / (1 + exp(-eta_MC))

# Accuracy scores (basis parameters)
acc_beta_BL <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_BL[j], var_BL[j])})
acc_beta_PG <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_PG[j], var_PG[j])})
acc_beta_PQ <- sapply(1:p, function(j) {accuracy(fit_MC$trace$beta[,j], mu_PQ[j], var_PQ[j])})
acc_beta_MC <- rep(1.0, times=p)

# Accuracy scores (linear predictor)
acc_eta_BL <- sapply(1:n, function(i) {accuracy(c(trace_eta_MC[,i]), eta_BL[i], var_eta_BL[i])})
acc_eta_PG <- sapply(1:n, function(i) {accuracy(c(trace_eta_MC[,i]), eta_PG[i], var_eta_PG[i])})
acc_eta_PQ <- sapply(1:n, function(i) {accuracy(c(trace_eta_MC[,i]), eta_PQ[i], var_eta_PQ[i])})
acc_eta_MC <- rep(1.0, times=n)

# Accuracy scores (predicted field)
loc <- expand.grid(x=lon_new, y=lat_new)
out <- sp::point.in.polygon(loc[,1], loc[,2], border[,1], border[,2])

acc_field_BL <- as.vector(as.matrix(psi_new %*% acc_beta_BL))
acc_field_PG <- as.vector(as.matrix(psi_new %*% acc_beta_PG))
acc_field_PQ <- as.vector(as.matrix(psi_new %*% acc_beta_PQ))
acc_field_MC <- NA

acc_field_BL[out == 0] <- NA
acc_field_PG[out == 0] <- NA
acc_field_PQ[out == 0] <- NA

# Posterior mean error (basis parameters)
err_mean_beta_BL <- abs(mu_BL - mu_MC)
err_mean_beta_PG <- abs(mu_PG - mu_MC)
err_mean_beta_PQ <- abs(mu_PQ - mu_MC)
err_mean_beta_MC <- rep(0.0, times=p)

# Posterior log-variance error (basis parameters)
err_var_beta_BL <- abs(log(var_BL) - log(var_MC))
err_var_beta_PG <- abs(log(var_PG) - log(var_MC))
err_var_beta_PQ <- abs(log(var_PQ) - log(var_MC))
err_var_beta_MC <- rep(0.0, times=p)

# Posterior mean error (linear predictor)
err_mean_eta_BL <- abs(eta_BL - eta_MC)
err_mean_eta_PG <- abs(eta_PG - eta_MC)
err_mean_eta_PQ <- abs(eta_PQ - eta_MC)
err_mean_eta_MC <- rep(0.0, times=n)

# Posterior standard deviation error (linear predictor)
err_var_eta_BL <- abs(log(var_eta_BL) - log(var_eta_MC))
err_var_eta_PG <- abs(log(var_eta_PG) - log(var_eta_MC))
err_var_eta_PQ <- abs(log(var_eta_PQ) - log(var_eta_MC))
err_var_eta_MC <- rep(0.0, times=n)

# Average accuracy/error scores
tvd_pdf_beta  <- rowMeans(rbind(1-acc_beta_BL, 1-acc_beta_PG, 1-acc_beta_PQ, 1-acc_beta_MC))
tvd_pdf_eta   <- rowMeans(rbind(1-acc_eta_BL, 1-acc_eta_PG, 1-acc_eta_PQ, 1-acc_eta_MC))
tvd_pdf_field <- rowMeans(rbind(1-acc_field_BL, 1-acc_field_PG, 1-acc_field_PQ, 1-acc_field_MC), na.rm=TRUE)
mae_mean_beta <- rowMeans(rbind(err_mean_beta_BL, err_mean_beta_PG, err_mean_beta_PQ, err_mean_beta_MC))
mae_var_beta  <- rowMeans(rbind(err_var_beta_BL, err_var_beta_PG, err_var_beta_PQ, err_var_beta_MC))
mae_mean_eta  <- rowMeans(rbind(err_mean_eta_BL, err_mean_eta_PG, err_mean_eta_PQ, err_mean_eta_MC))
mae_var_eta   <- rowMeans(rbind(err_var_eta_BL, err_var_eta_PG, err_var_eta_PQ, err_var_eta_MC))

# Summary table
summary <- data.frame(
  method = c("BL", "PG", "PQ", "MC"),
  niter = c(niter_BL, niter_PG, niter_PQ, niter_MC),
  exetime = round(c(exetime_BL, exetime_PG, exetime_PQ, exetime_MC), 3),
  timegain = round(1-c(exetime_BL, exetime_PG, exetime_PQ, exetime_MC)/exetime_MC, 4),
  elbo = round(c(elbo_BL, elbo_PG, elbo_PQ, elbo_MC), 3),
  tv_dist_field = round(tvd_pdf_field, 4),
  tv_dist_beta = round(tvd_pdf_beta, 4),
  mae_mean_beta = round(mae_mean_beta, 4),
  mae_var_beta = round(mae_var_beta, 4),
  tv_dist_eta = round(tvd_pdf_eta, 4),
  mae_mean_eta = round(mae_mean_eta, 4),
  mae_var_eta = round(mae_var_eta, 4),
  row.names = c(1:4)
)

print(summary)

if (SAVE) {
  filename <- "portland_mfvb_summary.csv"
  filepath <- paste(SAVEPATH, filename, sep="/")
  write.csv2(summary, file=filepath, row.names=FALSE)
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
  with(summary[1:3,], {
    barplt <- barplot(niter, names.arg=method, col=COLORS, border="white", xlab="", ylab="")
    text(barplt, niter-0.05*max(niter), niter)
    title(ylab="Iteration", main="Number of Iterations")
  })
  # exe-time
  with(summary[1:3,], {
    barplt <- barplot(exetime, names.arg=method, col=COLORS, border="white", xlab="", ylab="")
    text(barplt, exetime-0.05*max(exetime), round(exetime, 2))
    title(ylab="Time (s)", main="Execution Time")
  })
  dev.off()
}


if (FALSE) {
  filename <- "portland_mfvb_error_beta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(1,3))
  # TV distance
  df <- cbind("BL-VB"=1-acc_beta_BL, "PG-VB"=1-acc_beta_PG, "PQ-VB"=1-acc_beta_PQ)
  boxplot(df, col=c(2,3,4), ylim=c(0,1), main=expression(d[TV](q[VB], p[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - mu
  df <- cbind("BL-VB"=err_mean_beta_BL, "PG-VB"=err_mean_beta_PG, "PQ-VB"=err_mean_beta_PQ)
  boxplot(df, col=c(2,3,4), main=expression(MAE(mu[VB], mu[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - log(sigma)
  df <- cbind("BL-VB"=err_var_beta_BL, "PG-VB"=err_var_beta_PG, "PQ-VB"=err_var_beta_PQ)
  boxplot(df, col=c(2,3,4), main=expression(MAE(log~sigma[VB], log~sigma[MC])))
  abline(h=.0, col=7, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

if (FALSE) {
  filename <- "portland_mfvb_error_eta.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 4; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow=c(1,3))
  # TV distance
  boxplot(cbind("BL-VB"=1-acc_eta_BL, "PG-VB"=1-acc_eta_PG, "PQ-VB"=1-acc_eta_PQ), 
          col=c(2,3,4), ylim=c(0,1), main=expression(d[TV](q[VB], p[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - mu
  boxplot(cbind("BL-VB"=err_mean_eta_BL, "PG-VB"=err_mean_eta_PG, "PQ-VB"=err_mean_eta_PQ), 
          col=c(2,3,4), main=expression(MAE(mu[VB], mu[MC])))
  abline(h=.0, col=7, lty=1)
  # MAE - log(sigma)
  boxplot(cbind("BL-VB"=err_var_eta_BL, "PG-VB"=err_var_eta_PG, "PQ-VB"=err_var_eta_PQ), 
          col=c(2,3,4), main=expression(MAE(log~sigma[VB], log~sigma[MC])))
  abline(h=.0, col=7, lty=1)
  par(mfrow=c(1,1))
  dev.off()
}

### Spatial maps ----

if (FALSE) {
  filename <- "portland_mfvb_map_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  with(as.data.frame(locs), {
    par(mfrow = c(2,2))
    xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
    boundfun <- function() plot3D::lines2D(x=nodes[,1], y=nodes[,2], col="black", add=TRUE)
    plot3D::scatter2D(x=x, y=y, colvar=pr_BL, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="BL-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_PG, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PG-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_PQ, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="PQ-VB"); boundfun()
    plot3D::scatter2D(x=x, y=y, colvar=pr_MC, xlim=xlim, ylim=ylim, clim=clim, pch=19, xlab="", ylab="", main="MCMC"); boundfun()
    par(mfrow=c(1,1))
  })
  dev.off()
}

### Spatial fields ----

if (FALSE) {
  filename <- "portland_mfvb_field_pr.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 9; width <- 10; zoom <- 1
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow = c(2,2))
  xlim <- range(nodes[,1]); ylim <- range(nodes[,2]); clim <- c(0,1)
  plot_field(lon_new, lat_new, plogis(matrix(psi_new %*% mu_BL, ngrid, ngrid)), border, clim=c(0,1), main="BL-VB")
  plot_field(lon_new, lat_new, plogis(matrix(psi_new %*% mu_PG, ngrid, ngrid)), border, clim=c(0,1), main="PG-VB")
  plot_field(lon_new, lat_new, plogis(matrix(psi_new %*% mu_PQ, ngrid, ngrid)), border, clim=c(0,1), main="PQ-VB")
  plot_field(lon_new, lat_new, plogis(matrix(psi_new %*% mu_MC, ngrid, ngrid)), border, clim=c(0,1), main="MCMC")
  par(mfrow=c(1,1))
  dev.off()
}

if (FALSE) {
  filename <- "portland_mfvb_field_tvd.pdf"
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.5
  pdf(file=filepath, height=zoom*height, width=zoom*width)
  par(mfrow = c(1,3))
  xlim <- range(nodes[,1]); ylim <- range(nodes[,2])
  clim <- range(cbind(1-acc_beta_BL, 1-acc_beta_PG, 1-acc_beta_PQ, 1-acc_beta_MC))
  field_tvd_BL <- matrix(psi_new %*% (1-acc_beta_BL), ngrid, ngrid)
  field_tvd_PG <- matrix(psi_new %*% (1-acc_beta_PG), ngrid, ngrid)
  field_tvd_PQ <- matrix(psi_new %*% (1-acc_beta_PQ), ngrid, ngrid)
  plot_field(lon_new, lat_new, field_tvd_BL, border, clim=clim, main="BL-VB")
  plot_field(lon_new, lat_new, field_tvd_PG, border, clim=clim, main="PG-VB")
  plot_field(lon_new, lat_new, field_tvd_PQ, border, clim=clim, main="PQ-VB")
  par(mfrow=c(1,1))
  dev.off()
}

## PAPER PLOT ----

if (SAVE) {
  palette <- "inferno"
  filename <- paste0("portland_mfvb_field_tvd.pdf")
  filepath <- paste(IMGPATH, filename, sep="/")
  height <- 3; width <- 10; zoom <- 1.25

  loc <- expand.grid(x=lon_new, y=lat_new)
  out <- sp::point.in.polygon(loc[,1], loc[,2], border[,1], border[,2])
  dat <- data.frame(x=locs[,1], y=locs[,2], z=obs)
  tvd <- 1-cbind("BL"=acc_beta_BL, "PG"=acc_beta_PG, "PQ"=acc_beta_PQ)
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
