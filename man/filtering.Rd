\name{filtering}
\alias{filtering}
\alias{dsfilter}

\title{Double-Sided Symmetric Linear Filter}

\description{
Double-sided symmetric linear filter.
}

\usage{
filtering(x, mod, 
  trend = list(ar=1, ma=1, sigma2=NULL), 
  transitory = list(ar=1, ma=1, sigma2=NULL), 
  seasonal = list(ar=1, ma=1, sigma2=NULL),
  irregular.sigma2 = NULL,
  extend = 16, drift = FALSE)
dsfilter(x, w, mod, extend = 16)
}

\arguments{
\item{x}{a univariate time series.}
\item{mod}{an object of class \code{Arima}. See \code{\link[stats]{arima}}.}
\item{trend}{a list containing the coefficients and variance of the ARIMA model related to the 
trend component.}
\item{transitory}{a list containing the coefficients and variance of the ARIMA model related to the 
transitory component.}
\item{seasonal}{a list containing the coefficients and variance of the ARIMA model related to the 
seasonal component.}
\item{irregular.sigma2}{numeric, variance of the irregular component. If \code{NULL}, the 
estimate of the irregular component is not computed.}
\item{extend}{integer; if greater than zero, the series is extended by means of forecasts 
and backcasts based on the fitted model \code{mod}.}
\item{drift}{logical, if \code{TRUE} the \code{intercept} in the fitted model \code{mod} or an 
external regressor named \code{"drift"} is treated as a deterministic linear trend.}
\item{w}{a vector of filter coefficients (one side).}
}

\details{
These functions perform the convolution of the time series and the double-sided 
symmetric filter. They perform:

\code{stats::filter(c(rep(0, n-1), x, rep(0, n-1)), }

\code{+ filter=c(rev(w[-1]), w), method="convolution", sides=1)} 

where \code{n} is \code{length(x)}.

The design of the filter in the ARIMA-model-based decomposition procedure 
relies on the following result. 
The minimum mean squared error estimator of the component
is given by the ACGF of the model:
%
\deqn{%
\theta(L)x_t = \phi_n(L)\theta_s(L)a_t \,,}{%
theta(L)x[t] = phi_n(L)theta_s(L)a[t],}
%
where \eqn{\theta(L)}{theta(L)} is the MA of the model fitted to the observed data, 
\eqn{\theta_s(L)}{theta_s(L)} is the MA of the component (signal) to be estimated 
and \eqn{\phi_n(L)}{phi_n(L)} is the product of the AR polynomials 
of the remaining components.
The estimate of the signal, \eqn{\hat{s}_t}{s[t]}, is obtained by means of a 
double-sided symmetrical filter where the weights, \eqn{w}{w}, are the 
theoretical autocovariances of the model above:
%
\deqn{%
\hat{s}_t = \sum_{i=-\infty}^\infty w_i x_{t-i} \,.}{%
s[t] = sum_[i=-infty]^infty w[i] x[t-i].}
}

\value{
\code{filtering} returns a list of class \code{tsdecFilter} containing 
the series extended with forecasts (if \code{extend > 0})
(based on the ARMA model given as input),
the weights of one side of the filter for each component 
and the corresponding estimate of the components.

\code{dsfilter} returns the filtered time series.
}

%FIXME give reference
%\references{}

\seealso{
\code{\link[tsdecomp]{ARIMAdec}}, \code{\link[stats]{filter}}.
}

%\examples{}

\keyword{ts}
