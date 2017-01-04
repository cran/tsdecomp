\name{compare.acf}
\alias{compare.acf}
\alias{plot.tsdecAcf}

\title{Compare ACF of Theoretical, Estimator and Empirical Component}

\description{
Compute the AutoCorrelation functions of the following elements: 
the theoretical ARMA model of each component, the estimator for each component, 
the filtered or estimated components.
}

\usage{
compare.acf(x, mod, lag.max = 12, ...)
\method{plot}{tsdecAcf}(x, component = c("trend", "transitory", "seasonal"), ci = 0.95, 
  ci.type = c("ma", "white"), ci.class = c("estimator", "theoretical", "empirical"), 
  plot = TRUE, ...)
}

\arguments{
\item{x}{for \code{compare.acf}, an object of class \code{\link[tsdecomp]{ARIMAdec}};
for \code{plot.tsdecAcf}, an object of class \code{tsdecAcf} returned by \code{compare.acf}.}
\item{mod}{the object of class \code{Arima} decomposed in \code{x}. See \code{\link[stats]{arima}}.}
%
%\item{ar}{a list containing the AR polynomials of the model obtained for each component. 
%The output from \code{\link{roots.allocation}} is a valid input.}
%\item{ma}{a list containing the MA polynomials and variance of the the disturbance term 
%in the model obtained for each component. 
%The output from \code{\link{canonical.decomposition}} is a valid input.}
%\item{comp}{a matrix of time series objects containing the estimates of each component.
%The output from \code{\link{filtering}} is a valid input.}
%\item{diffpoly}{a named list containing the differencing polynomials to render the components 
%\code{comp} stationary. It is the element \code{ar$polys.nonstat} returned by \code{\link{roots.allocation}}.}
%
\item{lag.max}{maximum lag at which to calculate the autocorrelations.}
\item{component}{a character, the label of the component for which the ACF is to be obtained.}
\item{ci}{coverage probability for confidence interval. If this is zero or negative, 
confidence intervals are not computed}
\item{ci.type}{a character, the type of confidence interval. See details.}
\item{ci.class}{a character, the element that is taken as reference to computed the 
confidence intervals. Ignored if \code{ci.class='white'}.}
\item{plot}{logical, if \code{TRUE} the ACF is plotted.}
\item{\dots}{further arguments to be passed to \code{\link[stats]{acf}} 
and \code{\link[graphics]{plot}}.}
}

\details{
The ACF is obtained upon the stationary transformation of the models 
for the components and the estimators; i.e.,
non-stationary roots (if any) are removed from the AR polynomials.
The estimated components are also transformed according to the 
polynomials \code{x$ar$polys.nonstationary} that render the signals stationary.

Argument \code{ci.type} behaves similarly to the same argument in \code{\link[stats]{plot.acf}}.
If \code{ci.type = "white"}, the confidence bands are fixed to 
\eqn{t_{\alpha/2}/\sqrt(n)}{t_{alpha/2}/sqrt(n)}, where \eqn{n} is the number of observations 
in the fitted model \code{model}.
If \code{ci.type = "ma"}, confidence bands are obtained upon Bartlett's approximations 
for the standard deviations of the autocorrelations.
}

\value{
\code{compare.acf} returns the ACF of the components, respectively 
for their theoretical ARMA model, estimator and estimates.

\code{plot.tsdecAcf} displays a plot and returns a \code{\link[base]{invisible}}
copy of a matrix containing the confidence intervals.
}

%\references{}

\seealso{
\code{\link[tsdecomp]{ARIMAdec}}.
}

\examples{
# Airlines model and monthly data
y <- log(AirPassengers)
fit <- arima(y, order=c(0,1,1), seasonal=list(order=c(0,1,1)))
dec <- ARIMAdec(y, fit, extend=72)
cacf <- compare.acf(x = dec, mod=fit, lag.max=24)
plot(cacf, component="seasonal")
# unexpected discrepancy between the ACF of the estimator and the 
# ACF of the empirical signal
plot(cacf, component="trend")

# Nile time series
y <- Nile
fit <- arima(y, order=c(0,1,1))
dec <- ARIMAdec(y, fit, extend=16)
cacf <- compare.acf(x = dec, mod=fit, lag.max=24)
plot(cacf, component="trend")
}

\keyword{ts}
