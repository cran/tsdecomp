\name{ARMAacov}
\alias{ARMAacov}

\title{Compute Theoretical Autocovariances of an ARMA Model}

\description{
Compute the theoretical autocovariances of an ARMA model.
}

\usage{
ARMAacov(ar = numeric(0), ma = numeric(0), lag.max = max(p, q + 1), 
  sigma2 = 1)
}

\arguments{
\item{ar}{numeric vector of AR coefficients.}
\item{ma}{numeric vector of MA coefficients.}
\item{lag.max}{integer, maximum lag to be computed. The default is \code{max(p, q+1)},
where \code{p} and \code{q} are orders of the AR and MA terms,
\code{length(ar)} and \code{length(ma)}, respectively.
}
\item{sigma2}{numeric, the variance of the innovations.}
}

%\details{}

\value{
A vector of autocovariances named by lag order.
}

\note{
Based on \code{\link[stats]{ARMAacf}}.
}

\references{
Brockwell, P. J. and Davis, R. A. (1991) 
\emph{Time Series: Theory and Methods}, 
Second Edition. Springer.
\doi{10.1007/978-1-4419-0320-4}

Pollock, D. S. G. (1999)
\emph{A Handbook of Time-Series Analysis Signal Processing and Dynamics}. 
Academic Press. Chapter 17.
\doi{10.1016/B978-012560990-6/50002-6}
}

\seealso{
\code{\link[stats]{ARMAtoMA}}.
}

\examples{
# Autocovariances of an ARMA(2,1)
# method 1: using ARMAacov()
a1 <- ARMAacov(ar=c(0.8,-0.6), ma=0.4, lag.max=10)

# method 2: upon the coefficients of the infinite MA representation
psi <- c(1, ARMAtoMA(ar=c(0.8,-0.6), ma=0.4, lag.max=50))
a2 <- c(sum(psi^2), rep(0, length(a1)-1))
for (i in seq_along(a2[-1]))
  a2[i+1] <- sum(psi[seq_len(length(psi)-i)] * psi[-seq_len(i)])

# for a high enough number of 'psi' coefficients
# both methods are equivalent
all.equal(a1, a2, check.names=FALSE)
#[1] TRUE
}

\keyword{ts}
