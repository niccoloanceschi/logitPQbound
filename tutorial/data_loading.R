
get_std <- function(x,  weights = rep(1, nrow(x))){
  
  # `glmnet` style
  
  weights <- weights/sum(weights)
  xm <- drop(t(weights) %*% x)
  
  xv <- drop(t(weights) %*% x^2) - xm^2
  xv[abs(xv) < 10 * .Machine$double.eps] <- 0
  
  xs <- sqrt(xv)
  
  const_vars <- xs == 0
  if(any(const_vars)){xs[const_vars] <- 1}
  
  return(list(mean=xm,sd=xs)) ## Then: x <- scale(x, xm, xs)
}

get_Simulated_corX <- function(n, p, imbalance_ratio = 0.1, rho = 0.3, block_size = 5, seed = 123,
                               rescale_data=T) {
  
  set.seed(seed)
  
  p = p + p%%block_size
  
  # Step 1: Create a sparse correlation matrix (block structure)
  Sigma <- matrix(rho, nrow = block_size, ncol = block_size)
  diag(Sigma) <- 1
  cholSigma = chol(Sigma) # Cholesky decomposition to introduce correlation
  
  # Step 2: Generate uncorrelated normal data
  X_uncorrelated <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Step 3: Apply correlation structure to generate correlated data
  X <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:(p%/%block_size)) {
    idx = (i-1)*block_size + 1:block_size
    X[,idx] <- X_uncorrelated[,idx] %*% cholSigma 
  }
  
  # Step 4: Simulate coefficients
  beta <- runif(p,-3,3)
  
  # Step 4.1: Add intercept
  beta <- c(0,beta)
  X <- cbind(rep(1,n),X)
  p <- p+1
  
  # Step 5: Compute linear predictor and probabilities
  eta0 <- X %*% beta / sqrt(p)
  eta <- eta0 + LaplacesDemon::logit(imbalance_ratio)
  
  # Step 6: Compute probabilities
  pi <- LaplacesDemon::invlogit(eta)
  
  # Step 7: Sample binary outcomes
  y <- rbinom(n, size = 1, prob = pi)
  
  if(rescale_data){
    meansd <- get_std(X[,-1])
    X      <- cbind(rep(1,n),scale(X[,-1], meansd$mean, meansd$sd))
  }
  
  output = list(X = X, y = y, beta0 = beta, p=p, n=n)
  
  rm(X,y,eta,beta,eta0,pi)
  
  return(output)
}

get_Simulated <- function(n,p,n_test,imbalance_ratio,seed=1,rescale_data=T){
  
  set.seed(seed)
  
  beta0 <- c(0,runif(p-1,-5,5))
  
  X <- cbind(1,matrix(runif(n*(p-1)), ncol=p-1))
  y <- rbinom(n,1,LaplacesDemon::invlogit(X%*%beta0/sqrt(p)+LaplacesDemon::logit(imbalance_ratio)))
  
  X_test <- cbind(1,matrix(runif(n_test*(p-1)), ncol=p-1))
  y_test <- rbinom(n,1,LaplacesDemon::invlogit(X_test%*%beta0/sqrt(p)+LaplacesDemon::logit(imbalance_ratio)))
  
  if(rescale_data){
    meansd <- get_std(X[,-1])
    X      <- cbind(rep(1,n),scale(X[,-1], meansd$mean, meansd$sd))
    X_test <- cbind(rep(1,n_test),scale(X_test[,-1], meansd$mean, meansd$sd))
  }
  
  output = list(X=X,y=y,X_test=X_test,y_test=y_test,beta0=beta0,p=p,n=n,n_test=n_test)
  
  rm(X,X_test,y,y_test,beta0)
  
  return(output)
}

get_Leukemia <- function(n_frac=0.7,rescale_data=T,seed=1){
  
  set.seed(seed)
  
  load("~/Documents/GitHub/logitEN/data/Leukemia.RData")
  
  X <- as.matrix(Leukemia$x)
  y <- Leukemia$y
  
  n = floor(n_frac*length(y))
  n_test = length(y) - n
  p <- ncol(X)+1
  
  idx_test <- sample(1:length(y),n_test,replace=F)
  
  X_test <- cbind(1,X[idx_test,])
  y_test <- y[idx_test]
  
  X <- cbind(1,X[-idx_test,])
  y <- y[-idx_test]
  
  if(rescale_data){
    meansd <- get_std(X[,-1])
    X      <- cbind(rep(1,n),scale(X[,-1], meansd$mean, meansd$sd))
    X_test <- cbind(rep(1,n_test),scale(X_test[,-1], meansd$mean, meansd$sd))
  }
  
  output = list(p=p,n=n,n_test=n_test,X=X,y=y,X_test=X_test,y_test=y_test)
  
  rm(Leukemia,X,X_test,y,y_test)
  
  return(output)
}

get_Alzheimer <- function(n_frac=0.2,rescale_data=T,seed=1){
  
  load("~/Documents/GitHub/logitEN/data/Alzheimer_Interactions.RData")
  
  # cnt_u <- apply(X, 2, function(x) length(unique(x)))
  # X <- cbind(1,X[,cnt_u>5])
  
  X_test = X[-trainingSet,]
  y_test = y[-trainingSet]
  
  X = X[trainingSet,]
  y = y[trainingSet]
  
  if(n_frac<1){
    set.seed(seed)
    subSample <- sample(1:nrow(X),round(n_frac*nrow(X)))
    
    X = X[subSample,]
    y = y[subSample]
  }
  
  p = ncol(X)
  n = nrow(X)
  n_test = dim(X_test)[1]
  
  if(rescale_data){
    meansd <- get_std(X[,-1])
    X      <- cbind(rep(1,n),scale(X[,-1], meansd$mean, meansd$sd))
    X_test <- cbind(rep(1,n_test),scale(X_test[,-1], meansd$mean, meansd$sd))
  }
  
  output = list(p=p,n=n,n_test=n_test,X=X,y=y,X_test=X_test,y_test=y_test)
  
  rm(trainingSet,X,X_test,y,y_test)
  
  return(output)
}

get_PRISM <- function(n_frac=0.7,rescale_data=T,seed=1){
  
  Data <- readRDS('/Users/nico/Desktop/StelzerEGA_Data_cleaning/PRISM.rds')
  
  X <- cbind(rep(1,Data$n),Data$X_m[[1]])
  for(m in 2:Data$M){
    X <- cbind(X,Data$X_m[[m]])
  }
  y <- Data$yTrain
  
  if(n_frac<1){
    set.seed(seed)
    subSample <- sample(1:nrow(X),round(n_frac*nrow(X)))
    
    X = X[subSample,]
    y = y[subSample]
  }
  
  p = ncol(X)
  n = nrow(X)
  
  if(rescale_data){
    meansd <- get_std(X[,-1])
    X      <- cbind(rep(1,n),scale(X[,-1], meansd$mean, meansd$sd))
  }
  
  output = list(p=p,n=n,X=X,y=y)
  
  rm(X,y)
  
  return(output)
}


































































