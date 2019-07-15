\name{VAR_MLE}
\alias{VAR_MLE}
\title{
  Likelihood Based Variance Estimation in High-Dimensional Linear Models
}
\description{
  Use the maximum likelihood based method of Dicker and Erdogdu (2016) to estimate the residual variance in high-dimensional linear model.
}

\usage{VAR_MLE(y,x,max.iter=50,tol=1e-4)
}

\arguments{
  \item{y}{A length \eqn{n} vector of response.}
  
  \item{x}{A \eqn{n\times p} numeric design matrix for the model.}
  
  \item{max.iter}{Maximum number of iterations. Default is 50.}
  
  \item{tol}{Convergence threshhold. Default is 1e-4.}
}

\value{
  \item{sigma2}{The estimation of the residual variance.}
}


\references{Dicker, L. H. and Erdogdu, M. A. (2016). Maximum likelihood for variance estimation in high-dimensional linear models. In Proceedings of the 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016), 159-167. JMLR Workshop \eqn{\&} Conference Proceedings.
}

\examples{ 
n   <- 60
p   <- 100
beta <- c(sqrt(0.1/p)*rep(1,p/2),rep(0,p/2))
eps <- rnorm(n)
x <- matrix(rnorm(n*p),n,p)
y <- x\%*\%beta+eps
fit <- VAR_MLE(y,x)
}