library(glmnet)
library(MASS)
library(randomForest)

L2_Test <- function(x, z, y, estimate.method, split.num=1, is.PE=FALSE, threshold.type=NULL, 
                    lambda=NULL, boots.num=NULL, a.np=NULL, seed.fix) {
  
  # the sample size and dimension
  N <- length(y)
  p1 <- ncol(x); p2 <- ncol(z); p <- p1 + p2
  
  # p-value based on single data-splitting
  Pval_Single <- function(seed) {
    set.seed(seed + seed.fix)
    
    # data-splitting
    n <- floor(N/2)
    ind1 <- sample(1:N, n, replace=F)
    ind2 <- (1:N)[-ind1]
    
    Tn_Std <- function(ind1, ind2) {
      x1 <- x[ind1,,drop=F]; x2 <- x[ind2,,drop=F]
      z1 <- z[ind1,]; z2 <- z[ind2,]
      y1 <- y[ind1]; y2 <- y[ind2]
      
      # the prediction of g(Z)
      if(estimate.method == "Lasso") {
        # Lasso by R-package glmnet
        cv.model <- cv.glmnet(z1, y1, nfolds=10, intercept=FALSE)
        g.pred <- predict(cv.model, z2, s="lambda.min")
      }
      
      if(estimate.method == "RF") {
        # Random Forest by R-package randomForest
        data1 <- data.frame(y=y1, z1)
        data2 <- data.frame(y=y2, z2)
        rf.model <- randomForest(y ~ ., data=data1)
        g.pred <- predict(rf.model, data2)
      }
      
      # construction of test statistics
      x.mat <- x2 %*% t(x2)
      tr.Sigx2.hat <- (sum(x.mat^2) - sum(diag(x.mat^2)))/(n*(n-1))
      
      err.mat <- outer(as.numeric(y2-g.pred), as.numeric(y2-g.pred), "*")
      Tn <- (sum(err.mat*x.mat) - sum(diag(err.mat*x.mat)))/n
      sigma2.hat <- mean((y2-g.pred)^2)
      Tn.std <- Tn/(sigma2.hat*sqrt(2*tr.Sigx2.hat))
      
      if(is.PE) {
        # power enhancement
        Tk_Marginal <- function(k=1) {
          cat(paste0("Generating the ", k, "th marginal statistic...", "\r"))
          xk.mat <- outer(as.numeric(x2[,k]), as.numeric(x2[,k]), "*")
          Tk <- (sum(err.mat*xk.mat) - sum(diag(err.mat*xk.mat)))/(n*(n-1))
          if(threshold.type == "hard") {
            return(Tk)
          }
          
          if(threshold.type == "soft") {
            Tk.boots <- vector(length=boots.num)
            for (i in 1:boots.num) {
              set.seed(i + seed.fix)
              e.norm <- rnorm(n)
              e.mat <- outer(e.norm, e.norm, "*")
              Tk.H0 <- (sum(err.mat*xk.mat*e.mat) - sum(diag(err.mat*xk.mat*e.mat)))/(n*(n-1))
              Tk.boots[i] <- Tk.H0
            }
            Tk.boots.max <- max(Tk.boots)
            return(c(Tk, Tk.boots.max))
          }
        }
        
        Tk.all.ls <- lapply(1:p1, function(a){ Tk_Marginal(k=a) })
        Tk.all <- Reduce("rbind", Tk.all.ls)
        
        if(threshold.type == "hard") {
          delta.hard <- lambda*log(log(n))*(log(p1))^2/n
          T0.hard <- a.np*sum(abs(Tk.all)*ifelse((abs(Tk.all)-delta.hard) > 0, 1, 0))
          Tn.PE <- Tn.std + T0.hard
        }
        
        if(threshold.type == "soft") {
          T0.soft <- a.np*sum(abs(Tk.all[,1])*ifelse((abs(Tk.all[,1])-max(Tk.all[,2])) > 0, 1, 0))
          Tn.PE <- Tn.std + T0.soft
        }
        return(Tn.PE)
      }
      
      if(!is.PE) {
        return (Tn.std)
      }
    }
    
    Tn1.std <- Tn_Std(ind1, ind2)
    Tn2.std <- Tn_Std(ind2, ind1)
    
    Tn.crossfit <- (Tn1.std+Tn2.std)/sqrt(2)
    pval <- 1 - pnorm(Tn.crossfit)
    
    return(pval)
  }
  
  # p-value based on multi data-splitting
  pvals.single <- vector(length=split.num)
  for (i in 1:split.num) {
    cat(paste0("Performing the ", i, "th split...", "\n"))
    pvals.single[i] <- Pval_Single(seed=i)
  }
  
  if(split.num == 1) {
    pval.multi <- pvals.single
  } else {
    # combination of p-values by Meinshausen et al. (2009)
    Multi_Split_Pval <- function(pval) {
      gamma <- seq(0.05, 1, by=0.001)
      Q <- c()
      for (i in 1:length(gamma)) {
        Q[i] <- min(1, quantile(pval/gamma[i], gamma[i]))
      }
      Q1 <- min(1, (1-log(0.05))*min(Q))
      return(Q1)
    }
    pval.multi <- Multi_Split_Pval(pvals.single)
  }
  
  return(pval.multi)
}


## Example
# generate data
Gene_Data <- function(N, p1, p2, s1, c1, s2, c2, seed.rand) {
  set.seed(seed.rand)
  p <- p1 + p2
  beta.true <- c(rep(c1, s1), rep(0, p1-s1))
  gamma.true <- c(rep(c2, s2), rep(0, p2-s2))
  Sig <- toeplitz(0.5^seq(0, p-1))
  
  v <- mvrnorm(N, mu=rep(0, p), Sigma=Sig)
  x <- v[,1:p1]; z <- v[,-(1:p2)]
  g <- z%*%gamma.true
  error <- rnorm(N)
  y <- x%*%beta.true + g + error
  
  return(list(x=x, z=z, y=y))
}

data.H0 <- Gene_Data(N=200, p1=500, p2=500, s1=100, c1=0, s2=20, c2=0.5, seed.rand=13579)
data.H1 <- Gene_Data(N=200, p1=500, p2=500, s1=1, c1=1, s2=20, c2=0.5, seed.rand=13579)

## test H0: beta=0
# single data-splitting procedure
pval.single.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=1, seed.fix=0218)
pval.PE.hard.single.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="hard", lambda=0.9, a.np=5, seed.fix=0218)
pval.PE.soft.single.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="soft", boots.num=30, a.np=5, seed.fix=0218)
# # multiple data-splitting procedure
# pval.multi.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=30, seed.fix=0218)
# pval.PE.hard.multi.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=30, 
#                                   is.PE=TRUE, threshold.type="hard", lambda=0.9, a.np=5, seed.fix=0218)
# pval.PE.soft.multi.H0 <- L2_Test(x=data.H0$x, z=data.H0$z, y=data.H0$y, estimate.method="Lasso", split.num=30, 
#                                   is.PE=TRUE, threshold.type="soft", boots.num=30, a.np=5, seed.fix=0218)

## test H1: beta≠0
# single data-splitting procedure
pval.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, seed.fix=0218)
pval.PE.hard.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="hard", lambda=0.9, a.np=5, seed.fix=0218)
pval.PE.soft.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="soft", boots.num=30, a.np=5, seed.fix=0218)
# # multiple data-splitting procedure
# pval.multi.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=30, seed.fix=0218)
# pval.PE.hard.multi.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=30, 
#                                   is.PE=TRUE, threshold.type="hard", lambda=0.9, a.np=5, seed.fix=0218)
# pval.PE.soft.multi.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=30, 
#                                   is.PE=TRUE, threshold.type="soft", boots.num=30, a.np=5, seed.fix=0218)


