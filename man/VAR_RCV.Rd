\name{VAR_RCV}
\alias{VAR_RCV}
\title{
  Refitted Cross-Validation Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the refitted cross-validation based method of Fan et al. (2012) to estimate the residual variance in high-dimensional linear model.
}

\usage{VAR_RCV(y,x)
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
}

\value{
  \item{sigma2}{The estimation of the residual variance.}
}


\references{Fan, J., Guo, S. and Hao, N. (2012). Variance estimation using refitted cross-validation in ultrahigh-dimensional regression. Journal of Royal Statistical Society, Series B {\bold{74}}, 37-65.
}

\examples{ 
n   <- 80
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
eps <- rnorm(n)
x <- matrix(rnorm(n*p),n,p)
y <- x\%*\%beta+eps
fit <- VAR_RCV(y,x)
}