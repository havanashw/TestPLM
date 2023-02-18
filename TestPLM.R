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
