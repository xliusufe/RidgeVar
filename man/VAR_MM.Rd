\name{VAR_MM}
\alias{VAR_MM}
\title{
  Moments Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the moments based method of Dicker (2014) to estimate the residual variance in high-dimensional linear model.
}

\usage{VAR_MM(y,x,identity=F,Sigma=NULL)
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{identity}{Logical indicating  that the covariance of X is identity if identity=TRUE, and not if identity=FALSE. Default is FALSE. It is invalid if Sigma is not NULL.}
  
  \item{Sigma}{A \eqn{p\times p} matrix, which is the covariance of X. If Sigma=NULL, the sample covarince is given. Default is NULL.}
}

\value{
  \item{sigma2}{The estimator of the residual variance.}
}


\references{Dicker, L. H. (2014). Variance estimation in high-dimensional linear models. Biometrika {\bold{101}}, 269-284.
}

\examples{ 
n   <- 80
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
eps <- rnorm(n)
x <- matrix(rnorm(n*p),n,p)
y <- x\%*\%beta+eps
fit <- VAR_MM(y,x)
}