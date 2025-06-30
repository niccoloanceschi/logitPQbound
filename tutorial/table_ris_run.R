
# Imports ----------------------------------------------------------------------

library(viridis)
library(microbenchmark)
library(glmnet)
library(logitPQbound)
library(dplyr)
library(ggplot2)
library(rstudioapi)

source(file.path(dirname(getSourceEditorContext()$path),'data_loading.R'))

options(digits=5)

# Hyper-Parameters -------------------------------------------------------------

folder_out <- file.path(dirname(getSourceEditorContext()$path),'ris_paper/')
if (!dir.exists(folder_out)) {
  dir.create(folder_out, recursive = TRUE)
}

maxIter = 1e4 # Maximum number of iterations
tol=1e-10 # Threshold at which the algorithm stops
maxIterIn = 1e2 # Maximum number of inner iterations for PQ bound optimization
tolIn=1e-1 # Threshold at which the inner optimization of PQ bound stops
eps=1e-7 # Hessian's diagonal extra term for stability

rescaled <- T

my_seed = 123 

nRep=2 # 10 # 1 # Number of replicate runs

which_lambda <- 'Custom' # in c('Custom', 'Scaled') 

if (which_lambda=='Custom'){
  log_lambda <- 0 #
  lambda  <- exp(log_lambda)
  l_print <- sprintf("%.2f",log_lambda)
} else if (which_lambda=='Scaled'){
  lFrac = 2
  l_print <- paste0('_p_25e',lFrac)
  # TO BE COMPUTED CASE BY CASE
  # lambda  <- ncol(Data$X)/(25*10^lFrac)
} else {
  stop("Select lambda via between 'Custom', 'Scaled'")
}

# Select dataset ---------------------------------------------------------------

## Import data

which_data = 'Leukemia' # in c('Leukemia', 'Alzheimer', 'Random', 'PRISM')

if (which_data == 'Leukemia'){
  Data <- get_Leukemia(n_frac=0.49,rescale_data=rescaled,seed=my_seed)
} else if (which_data == 'Alzheimer'){
  Data <- get_Alzheimer(n_frac=0.15,rescale_data=rescaled,seed=my_seed)
} else if (which_data == 'PRISM'){
  Data <- get_PRISM(n_frac=0.37,rescale_data=rescaled,seed=my_seed)
} else if (which_data == 'Random') {
  small_n_large_p=T
  if(small_n_large_p){
    n=50
    p=200
  } else {
    n=200
    p=50
  }
  nTest=50
  r=0.7
  # Data <- get_Simulated(n,p,nTest,imbalance_ratio=r,seed=my_seed)
  Data <- get_Simulated_corX(n=n,p=p,imbalance_ratio=r,seed = my_seed)
} else {
  stop("Choose data set between 'Leukemia' and 'Alzheimer'")
}

beta_init = rep(0,Data$p)

if(which_lambda == 'Scaled') {
  lambda <- Data$p/(25*10^lFrac)
} else if (which_lambda == 'CV') {
  lambda_cv_glmnet <- cv.glmnet(Data$X[,-1],Data$y,alpha=alpha,family="binomial",
                                standardize=F, thres=tol, maxit=maxit)$lambda.min
  lambda <- round(n*lambda_cv_glmnet,2)
}

dim(Data$X)

# Blocked Optimization [ridge only] --------------------------------------------

## All lambda ------------------------------------------------------------------

if(F){
  lambda_all <- exp(seq(10,-5,length.out=100))
  lambda_all <- sort(lambda_all,decreasing=F)
  
  loglik_NR <- rep(NA,length(lambda_all))
  for(l in 1:length(lambda_all)){
    fit_NR_bl  <- ridge_logit(Data$X, Data$y, type='NR', lambda=lambda_all[l], beta_start=beta_init, 
                              maxiter=maxIter, tol=tol)
    loglik_NR[l] <- fit_NR_bl$Convergence[nrow(fit_NR_bl$Convergence),2]
  }
  graphics.off()
  plot(log(lambda_all),loglik_NR, type='l',ylab='loglik',xlab='log(lambda)',col='#1B420D')
}

## Evaluation ------------------------------------------------------------------

time_rep <- matrix(0,nRep+1,5)

for(rr in 1:(nRep+1)){
  # ----- NR ----- #
  ptm <- proc.time()
  fit_NR <- ridge_logit(Data$X, Data$y, type='NR', lambda=lambda, beta_start=beta_init, 
                           maxiter=maxIter, tol=tol)
  end_time <- proc.time() - ptm
  time_rep[rr,1] <- end_time["user.self"] + end_time["sys.self"]
  
  # ----- BL ----- #
  ptm <- proc.time()
  fit_BL <- ridge_logit(Data$X, Data$y, type='BL', lambda=lambda, beta_start=beta_init, 
                           maxiter=maxIter, tol=tol)
  end_time <- proc.time() - ptm
  time_rep[rr,2] <- end_time["user.self"] + end_time["sys.self"]
  
  # ----- PG ----- #
  ptm <- proc.time()
  fit_PG <- ridge_logit(Data$X, Data$y, type='PG', lambda=lambda, beta_start=beta_init, 
                           maxiter=maxIter, tol=tol)
  end_time <- proc.time() - ptm
  time_rep[rr,3] <- end_time["user.self"] + end_time["sys.self"]
  
  # ----- PQ ----- #
  ptm <- proc.time()
  fit_PQ <- ridge_logit(Data$X, Data$y, type='PQ', lambda=lambda, beta_start=beta_init, 
                            maxiter=maxIter, tol=tol, maxiter_in=maxIterIn, tol_in=tolIn)
  end_time <- proc.time() - ptm
  time_rep[rr,4] <- end_time["user.self"] + end_time["sys.self"]
  
  # ----- PQ-BOOST ----- #
  beta_init_plq <- beta_init
  beta_init_plq[1] <- 10*(1-2*(sum(Data$y)/length(Data$y)<0.5))
  
  ptm <- proc.time()
  fit_PQ_boost <- ridge_logit(Data$X, Data$y, type='PQ', lambda=lambda, beta_start=beta_init_plq, 
                                 maxiter=maxIter, tol=tol, maxiter_in=maxIterIn, tol_in=tolIn)
  end_time <- proc.time() - ptm
  time_rep[rr,5] <- end_time["user.self"] + end_time["sys.self"]
}


## Create Tables ---------------------------------------------------------------

if(nRep>1){
  time_avg <- apply(time_rep[-1,], 2, median)
} else {
  time_avg <- c(time_rep[2,])
}

iter_conv <- c(nrow(fit_NR$Convergence),nrow(fit_BL$Convergence),
                  nrow(fit_PG$Convergence),nrow(fit_PQ$Convergence),
                  nrow(fit_PQ_boost$Convergence))

ris_table <- data.frame(mthd=c('NR','BL','PG','PQ','PQ_boost'),
                           nIter=iter_conv , time=round(time_avg,2))

ris_full <- list('table'=ris_table, 'NR'=fit_NR, 'BL'=fit_BL,
               'PG'=fit_PG, 'PQ'=fit_PQ,
               'PQ_boost'=fit_PQ_boost)

ris_table

## Saving Outputs --------------------------------------------------------------

write.csv(ris_table, paste0(folder_out,"ris_",which_data,"_l",l_print,".csv"))

saveRDS(ris_full, file = paste0(folder_out,"ris_",which_data,"_l",l_print,".rds"))


















































        