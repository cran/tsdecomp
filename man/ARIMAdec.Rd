\name{ARIMAdec}
\alias{ARIMAdec}
\alias{plot.ARIMAdec}
\alias{print.ARIMAdec}

\title{ARIMA-Model-Based Decomposition of Time Series}

\description{
This is the main function for the ARIMA-model-based decomposition of a time series.
}

\usage{
ARIMAdec(x, mod, width = c(0.035, 0.035), min.modulus = 0.4, 
  extend = 16, drift = FALSE, optim.tol = 1e-04, ...)
\method{print}{ARIMAdec}(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
\method{plot}{ARIMAdec}(x, ...)
}

\arguments{
\item{x}{for \code{ARIMAdec}, a univariate time series; 
for \code{plot.ARIMAdec} and \code{print.ARIMAdec}, 
an object of class \code{ARIMAdec} returned by \code{ARIMAdec}.}
\item{mod}{an object of class \code{Arima}. See \code{\link[stats]{arima}}.}
%
\item{width}{numeric of length two, width of the interval of frequencies allocated to the trend 
and the seasonal components (measured in radians). If a numeric of length one is passed as argument, the same 
width is used for both components. See \code{\link{roots.allocation}}.}
\item{min.modulus}{numeric, minimum modulus of the roots assigned to the trend component. 
See \code{\link{roots.allocation}}.}
%
\item{extend}{integer; if greater than zero, the series is extended by means of forecasts 
and backcasts based on the fitted model \code{mod}.
See \code{\link{filtering}}.}
\item{drift}{logical, if \code{TRUE} the \code{intercept} in the fitted model \code{mod} or an 
external regressor named \code{"drift"} is treated as a deterministic linear trend.
See \code{\link{filtering}}.}
%
\item{optim.tol}{numeric, the convergence tolerance to be used by \code{\link[stats]{optimize}}.}
%
\item{units}{character, the units in which the argument of the roots are printed. \code{units="pi"} prints the 
argument in radians as multiples of \eqn{\pi}{pi}.}
\item{digits}{numeric, the number of significant digits to be used by \code{\link[base]{print}}.}
%
\item{\dots}{further arguments to be passed to \code{\link{poly2acgf}} or to 
\code{\link{plot.tsdecFilter}} and \code{\link[base]{print}} methods.}
}

\details{
This function is a wrapper to the sequence of calls to 
\code{\link{roots.allocation}}, \code{\link{pseudo.spectrum}}, 
\code{\link{canonical.decomposition}} and \code{\link{filtering}}.
}

\value{
An object of class \code{ARIMAdec} containing the following:
1) \code{ar}: the output from \code{{roots.allocation}},
2) \code{spectrum}: the output from \code{{pseudo.spectrum}},
3) \code{ma}:  the output from \code{{canonical.decomposition}},
4) \code{xextended}: the series extended with backcasts and forecasts (if \code{extend > 0}),
5) \code{filters}: the filters returned by \code{{filtering}},
6) \code{components}: the estimated components returned by \code{{filtering}}.
}

\references{
Burman, J. P. (1980) 
\sQuote{Seasonal Adjustment by Signal Extraction}.
\emph{Journal of the Royal Statistical Society. Series A (General)},
\bold{143}(3), pp. 321-337.
\doi{10.2307/2982132}

Hillmer, S. C. and Tiao, G. C. (1982) 
\sQuote{An ARIMA-Model-Based Approach to Seasonal Adjustment}.
\emph{Journal of the American Statistical Association},
\bold{77}(377), pp. 63-70.
\doi{10.1080/01621459.1982.10477767}
}

\seealso{
\code{\link{canonical.decomposition}},
\code{\link{filtering}},
\code{\link{pseudo.spectrum}},
\code{\link{roots.allocation}}.
}

\examples{
# Airlines model and monthly data
y <- log(AirPassengers)
fit <- arima(y, order=c(0,1,1), seasonal=list(order=c(0,1,1)))
dec <- ARIMAdec(y, fit, extend=72)
dec
plot(dec)

# JohnsonJohnson quarterly data
y <- log(JohnsonJohnson)
fit <- arima(y, order=c(0,1,1), seasonal=list(order=c(0,1,1)))
dec <- ARIMAdec(y, fit, extend=16)
dec
plot(dec)

# Nile annual data
# this series is better modelled as a level shift at 
# observation 29 and a mean (no ARMA structure),
# here the shift is ignored for illustration of the 
# decomposition of the fitted ARIMA(0,1,1) model
y <- Nile
fit <- arima(y, order=c(0,1,1))
dec <- ARIMAdec(y, fit, extend=72)
dec
plot(dec, overlap.trend=TRUE, args.trend=list(col="blue"))
}

\keyword{ts}
