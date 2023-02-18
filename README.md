# TestPLM
R code for "Tests for ultrahigh-dimensional partially linear regression models."

## Usage

```{R}
L2_Test(x, z, y, estimate.method, split.num = 1, is.PE = FALSE, threshold.type = NULL, 
        lambda = NULL, boots.num = NULL, a.np = NULL, seed.fix)
```

## Required Packages
- `glmnet`
- `MASS`
- `randomForest`

## Inputs:

- `x`: A matrix of n*p1, where n is the sample size and p1 is the dimension of interested covariates.
- `z`: A matrix of n*p2, where p2 is the dimension of potential nuisance variables.
- `y`: The response variable, with length n.
- `estimate.method`: The method used to estimate the high-dimensional nuisance function `g`. The options are `"Lasso"` or `"RF"`.
- `split.num`: The number of data-splitting. The default value is 1, indicating single data-splitting procedure is used.
- `is.PE`: A logical value specifying whether the power enhancement is used. The default value is `FALSE`.
- `threshold.type`: The type of thresholding to use in the power enhancement procedure. The options are `hard` or `soft`. If `hard`, the hard thresholding is used, and if `soft`, the soft thresholding is used.
- `lambda`: The tuning parameter for the power enhancement procedure with the hard threshold.
- `boots.num`: The tuning parameter for the power enhancement procedure with the soft threshold.
- `a.np`: The tuning parameter for the power enhancement procedure.
- `seed.fix`: The fixed seed.

## Examples:

```{R}
library(glmnet)
library(MASS)
library(randomForest)
source("TestPLM.R")

## generate data
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
                                  
## test H1: beta≠0
# single data-splitting procedure
pval.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, seed.fix=0218)
pval.PE.hard.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="hard", lambda=0.9, a.np=5, seed.fix=0218)
pval.PE.soft.single.H1 <- L2_Test(x=data.H1$x, z=data.H1$z, y=data.H1$y, estimate.method="Lasso", split.num=1, 
                                  is.PE=TRUE, threshold.type="soft", boots.num=30, a.np=5, seed.fix=0218)
```
