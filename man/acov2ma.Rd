\name{acov2ma}
\alias{acov2ma}
\alias{acov2ma.init}

\title{Convert Autocovariances to Coefficients of a Moving Average}

\description{
Convert autocovariances to coefficients of a moving average.
}

\usage{
acov2ma.init(x, tol = 0.00001, maxiter = 100)
acov2ma(x, tol = 1e-16, maxiter = 100, init = NULL)
}

\arguments{
\item{x}{a numeric vector containing the autocovariances.}
\item{tol}{numeric, convergence tolerance.}
\item{maxiter}{numeric, maximum number of iterations.}
\item{init}{numeric, vector of initial coefficients.}
}

\details{
\code{acov2ma.init} is based on procedure (17.35) described in Pollock (1999).
\code{acov2ma} is the Newton-Raphson procedure (17.39) 
described in the same reference.
}

\value{
A list containing the vector of coefficients and the variance of the 
innovations in the moving average model; 
convergence code and number of iterations.
}

\references{
Pollock, D. S. G. (1999)
\emph{A Handbook of Time-Series Analysis Signal Processing and Dynamics}. 
Academic Press. Chapter 17.
\doi{10.1016/B978-012560990-6/50002-6}
}

%\seealso{}

\examples{
set.seed(123)
x <- arima.sim(n=200, model=list(ma=c(0.7,-0.3)))
#sample autocovariances
a <- c(var(x), cov(x[-1], x[-200]), cov(x[-c(1,2)], x[-c(199,200)]))
#inferred coefficients and variance
acov2ma(a, init=acov2ma.init(a, maxit=10)$macoefs)
#compare with maximum-likelihood
arima(x, order=c(2,0,0), include.mean=FALSE)
}

\keyword{ts}
