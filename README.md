# RidgeVar
R package "RidgeVar" for estimation of error variance via ridge regression. Provide several methods to estimate the error variance for high-dimensional linear regression models, which includes the ridge regression based method of Liu et al. (2020), the refitted cross validation of Fan et al. (2012), the maximum likelihood based method of Dicker and Erdogdu (2016), and the moments based method of Dicker (2014).

# Installation

    #install.packages("devtools")
    library(devtools)
    install_github("xliusufe/RidgeVar")

# Usage

   - [x] [RidgeVar-manual.pdf](https://github.com/xliusufe/RidgeVar/blob/master/inst/RidgeVar-manual.pdf) ---------- Details of the usage of the package.
   - [x] [ridgevar](https://github.com/xliusufe/RidgeVarpy) ------------------------ the corresponding `Python` package
# Example
    library(RidgeVar)

    # Example 1. the usage of the function "VAR_RR" when `q=1`   
    n   <- 80
    p   <- 100
    beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
    eps <- rnorm(n)
    x <- matrix(rnorm(n*p),n,p)
    y <- x%*%beta+eps
    fit <- VAR_RR(y,x)

    Example 2. the usage of the function "VAR_RR" when `q>1`

    n   <- 80
    p   <- 100
    q   <- 5
    rho <- 0.5
    beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
    sig  <- chol( outer(1:q,1:q,function(i,j) rho^abs(i-j)) )
    y <- matrix(rnorm(n*q),n,q) %*% sig
    x <- matrix(rnorm(n*p),n,p)
    y[,1] <- y[,1] + x %*% beta
    fit <- VAR_RR(y,x)

    # Example 3. the usage of the function "COV_RR"
    
    n   <- 80
    p   <- 100
    q   <- 5
    rho <- 0.5
    beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
    sig  <- chol( outer(1:q,1:q,function(i,j) rho^abs(i-j)) )
    y <- matrix(rnorm(n*q),n,q) %*% sig
    x <- matrix(rnorm(n*p),n,p)
    y[,1] <- y[,1] + x %*% beta
    fit <- COV_RR(y,x)

    
# References
Dicker, L. H. (2014). Variance estimation in high-dimensional linear models.  Biometrika 101, 269-284.

Dicker, L. H. and Erdogdu, M. A. (2016). Maximum likelihood for variance estimation in high-dimensional linear models. In  Proceedings     of the 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016), 159-167. JMLR Workshop & Conference  Proceedings.

Fan, J., Guo, S. and Hao, N. (2012). Variance estimation using refitted cross-validation in ultrahigh-dimensional regression. Journal of Royal Statistical Society Series B 74, 37-65.

Liu, X., Zheng, S. and Feng, X. (2020). Estimation of error variance via ridge regression. Biometrika, 107, 481-488. [DOI: 10.1093/biomet/asz074](DOI:10.1093/biomet/asz074)

# Development
This R package is developed by Xu Liu (liu.xu@sufe.edu.cn) and Xiao Zhang.
