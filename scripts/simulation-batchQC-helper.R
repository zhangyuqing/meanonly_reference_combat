genewise_stats <- function(data, batch){
  nbatch <- length(unique(batch))
  stats.data <- list()
  library(moments)
  func_list <- c(mean, var, skewness, kurtosis)
  for(k in 1:length(func_list)){
    stats.data[[k]] <- matrix(0, nrow=nrow(data), ncol=nbatch,
                              dimnames=list(rownames(data), 
                                            paste("Batch", unique(batch), sep="")))
    for(i in 1:nbatch){
      stats.data[[k]][, paste("Batch", unique(batch)[i], sep="")] <- apply(data[, batch==unique(batch)[i]], 1, func_list[[k]])
    }
  }
  names(stats.data) <- c("Mean", "Variance", "Skew", "Kurtosis")
  return(stats.data)
}


#### batchQC_condition_adjusted from BatchQC utils.R under R/
batchQC_condition_adjusted <- function(data.matrix, batch, condition) {
  nlc <- nlevels(as.factor(condition))
  if (nlc <= 1)  {
    return(data.matrix)
  }
  P <- nlevels(as.factor(batch))
  if (P <= 1)  {
    pdata <- data.frame(condition)
    X <- model.matrix(~as.factor(condition), data = pdata)
    P <- 1
  } else  {
    pdata <- data.frame(batch, condition)
    X <- model.matrix(~as.factor(batch) + as.factor(condition), data=pdata)
  }
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(data.matrix))
  condition_adjusted <- data.matrix - t(as.matrix(X[, -c(1:P)]) %*% 
                                          beta[-c(1:P), ])
  return(condition_adjusted)
}


fitMoments <- function(x) {
  # Compute Mean and Variance
  smean <- apply(x, 2, mean)
  svariance <- apply(x, 2, var)
  sskew <- apply(x, 2, skewness)
  skurt <- apply(x, 2, kurtosis)
  sfit <- cbind(smean, svariance, sskew, skurt)
  colnames(sfit) = c("Mean", "Variance", "Skew", "Kurtosis")
  rownames(sfit) = colnames(x)
  return(sfit)
}


delta_f.pvalue <- function(dat, mod, mod0) {
  ## F-test (full/reduced model) and returns R2 values
  ## (full/reduced) as well.
  mod00 <- matrix(rep(1, ncol(dat)), ncol = 1)
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  
  resid <- dat - dat %*% mod %*% solve(t(mod) %*% mod) %*% t(mod)
  rss1 <- rowSums(resid * resid)
  rm(resid)
  
  if (df0 > 0)  {
    resid0 <- dat - dat %*% mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0)
  } else {
    resid0 <- dat
  }
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
  
  resid00 <- dat - dat %*% mod00 %*% solve(t(mod00) %*% mod00) %*% t(mod00)
  rss00 <- rowSums(resid00 * resid00)
  rm(resid00)
  
  r2_full <- 1 - rss1/rss00
  r2_reduced <- 1 - rss0/rss00
  
  #delta <- (apply(dat, 1, mean) * 0.01)^2
  delta <- n * (apply(dat, 1, mean) * 0.05)^2
  #delta <- 0
  p <- 1
  if (df1 > df0)  {
    fstats <- ((rss0 - rss1)/(df1 - df0))/(delta + rss1/(n - df1))
    p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  }
  return(list(p = p, r2_full = r2_full, r2_reduced = r2_reduced))
}
