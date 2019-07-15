\name{VAR_RR}
\alias{VAR_RR}
\title{
  Ridge Regression Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the ridge regression based method of Xu Liu et al. (2019) to estimate the residual variance in linear model. The proposal is valid under both low- and high-dimensional models, and performs well not only for the non-sparse cases but also for the sparse ones. 
}

\usage{VAR_RR(y,x,eta=NULL,alpha=0.1)
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{eta}{The tunning parameter using in ridge regression. Default is NULL.}

  \item{alpha}{A constant to justify the tunning parameter when eta=NULL. Default is 0.1.}
  
}

\value{
  \item{sigma2}{The estimation of the residual variance.}
}


\references{Liu, X., Zheng, S. and Feng, X. (2019). Estimation of error variance via ridge regression. Manuscript.
}

\examples{ 
n   <- 60
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
eps <- rnorm(n)
x <- matrix(rnorm(n*p),n,p)
y <- x\%*\%beta+eps
fit <- VAR_RR(y,x)
}