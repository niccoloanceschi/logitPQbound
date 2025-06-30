
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

folder_in <- file.path(dirname(getSourceEditorContext()$path),'ris_paper/')
folder_out <- paste0(folder_in,'plots/')
if (!dir.exists(folder_out)) {
  dir.create(folder_out, recursive = TRUE)
}

rescaled <- T

my_seed = 123 

which_lambda <- 'Custom' # in c('Custom', 'Scaled') 

if (which_lambda=='Custom'){
  log_lambda <- 0
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

## Loading Results -------------------------------------------------------------

ris_full = readRDS(file=paste0(folder_in,"ris_",which_data,"_l",l_print,".rds"))

## Plot Convergence ------------------------------------------------------------

file_out <- paste0(folder_out,"ris_",which_data,"_l",l_print,".pdf")
pdf(file=file_out,height=4,width=6)
par(mar=c(3.1, 5.2, 1.5, 1.5))

y_min = min(fit_BL$Convergence[-1,2])
y_max = max(fit_BL$Convergence[-1,2])

xMax <- 50 # 200

xMax_NR <- min(xMax,nrow(fit_NR$Convergence))
xMax_PQ <- min(xMax,nrow(fit_PQ_boost$Convergence))
xMax_PQ <- min(xMax,nrow(fit_PQ$Convergence))
xMax_PG <- min(xMax,nrow(fit_PG$Convergence))
xMax_BL <- min(xMax,nrow(fit_BL$Convergence))

p <- plot(c(2:xMax_PQ), fit_PQ_boost$Convergence[2:xMax_PQ,2], col='#B9DF92', pch=5, lwd=1, cex=0.5,
             ylab="Penalized Log-Likelihood", xlab='Iterations', xlim=c(1,xMax), ylim=c(y_min,y_max))
p <- p + points(c(2:xMax_PQ), fit_PQ$Convergence[2:xMax_PQ,2], col='#FFB000', pch=1, lwd=1, cex=0.5)
p <- p + points(c(2:xMax_PG), fit_PG$Convergence[2:xMax_PG,2], col='#DC267F', pch=2, lwd=1, cex=0.5)
p <- p + points(c(2:xMax_BL), fit_BL$Convergence[2:xMax_BL,2], col='#785EF0', pch=0, lwd=1, cex=0.5)
legend('bottomright', legend=c(expression("PQ"["BOOST"]), "PQ", "PG", "BL"),
       col=c("#B9DF92", "#FFB000", "#DC267F", "#785EF0"),
       pch=c(5,1,2,0), lwd=1.2, cex=0.8, lty=c(0,0,0))

dev.off()

## Adding Newton Raphson -------------------------------------------------------

file_out <- paste0(folder_out,"ris_",which_data,"_l",l_print,"_NR",".pdf")
pdf(file=file_out,height=4,width=6)
par(mar=c(3.1, 5.2, 1.5, 1.5))

y_min = min(fit_BL$Convergence[-1,2])
y_max = max(fit_BL$Convergence[-1,2])

p <- plot(c(2:xMax_NR), fit_NR$Convergence[2:xMax_NR,2], col='#1B420D', pch=6, lwd=1, cex=0.5,
             ylab="Penalized Log-Likelihood", xlab='Iterations', xlim=c(1,xMax), ylim=c(y_min,y_max))
p <- p + points(c(2:xMax_PQ), fit_PQ_boost$Convergence[2:xMax_PQ,2], col='#B9DF92', pch=5, lwd=1, cex=0.5)
p <- p + points(c(2:xMax_PQ), fit_PQ$Convergence[2:xMax_PQ,2], col='#FFB000', pch=1, lwd=1, cex=0.5)
p <- p + points(c(2:xMax_PG), fit_PG$Convergence[2:xMax_PG,2], col='#DC267F', pch=2, lwd=1, cex=0.5)
p <- p + points(c(2:xMax_BL), fit_BL$Convergence[2:xMax_BL,2], col='#785EF0', pch=0, lwd=1, cex=0.5)
legend('bottomright', legend=c("NR", expression("PQ"["BOOST"]), "PQ", "PG", "BL"),
       col=c("#1B420D", "#B9DF92", "#FFB000", "#DC267F", "#785EF0"),
       pch=c(6,5,1,2,0), lwd=1.2, cex=0.8, lty=c(0,0,0))

dev.off()

## Log Scale Diff -----


## Plot Convergence ------------------------------------------------------------

file_out <- paste0(folder_out,"ris_",which_data,"_l",l_print,"_log_diff.pdf")
pdf(file=file_out,height=4,width=6)
par(mar=c(3.1, 5.2, 1.5, 1.5))

xMax <- max(nrow(fit_PQ_boost$Convergence),nrow(fit_PQ$Convergence),
            nrow(fit_PG$Convergence),nrow(fit_BL$Convergence))

xMax_PQ_boost <- min(xMax,nrow(fit_PQ_boost$Convergence))
xMax_PQ <- min(xMax,nrow(fit_PQ$Convergence))
xMax_PG <- min(xMax,nrow(fit_PG$Convergence))
xMax_BL <- min(xMax,nrow(fit_BL$Convergence))

y_PQ_boost = - fit_PQ_boost$Convergence[2:(xMax_PQ_boost-1),2] + fit_PQ_boost$Convergence[3:xMax_PQ_boost,2]
y_PQ = - fit_PQ$Convergence[2:(xMax_PQ-1),2] + fit_PQ$Convergence[3:xMax_PQ,2]
y_PG = - fit_PG$Convergence[2:(xMax_PG-1),2] + fit_PG$Convergence[3:xMax_PG,2]
y_BL = - fit_BL$Convergence[2:(xMax_BL-1),2] + fit_BL$Convergence[3:xMax_BL,2]

y_min = min(abs(c(y_PQ_boost,y_PQ,y_PG,y_BL)))
y_max = max(y_PQ_boost,y_PQ,y_PG,y_BL)

p <- plot(c(2:(xMax_PQ_boost-1)), y_PQ_boost, col='#B9DF92', pch=5, lwd=1, cex=0.5, log='y',
          ylab="Sequential Differences in Target", xlab='Iterations', xlim=c(1,xMax), ylim=c(y_min,y_max))
p <- p + points(c(2:(xMax_PQ-1)), y_PQ, col='#FFB000', pch=1, lwd=1, cex=0.5)
p <- p + points(c(2:(xMax_PG-1)), y_PG, col='#DC267F', pch=2, lwd=1, cex=0.5)
p <- p + points(c(2:(xMax_BL-1)), y_BL, col='#785EF0', pch=0, lwd=1, cex=0.5)
legend('topright', legend=c(expression("PQ"["BOOST"]), "PQ", "PG", "BL"),
       col=c("#B9DF92", "#FFB000", "#DC267F", "#785EF0"),
       pch=c(5,1,2,0), lwd=1.2, cex=0.8, lty=c(0,0,0))

dev.off()


