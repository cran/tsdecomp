\name{roots.allocation}
\alias{roots.allocation}
\alias{plot.tsdecARroots}
\alias{print.tsdecARroots}

\title{Allocation of Autoregressive Roots}

\description{
Allocate the roots of the autoregressive polynomial from a fitted ARIMA model to 
trend, transitory and seasonal components.
}

\usage{
roots.allocation(x, width = c(0.035, 0.035), min.modulus = 0.4)
\method{plot}{tsdecARroots}(x, xlim, ylim, ...)
\method{print}{tsdecARroots}(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
}

\arguments{
\item{x}{for \code{roots.allocation}, an object of class \code{Arima} (see \code{\link[stats]{arima}});
for \code{print.tsdecARroots} and \code{plot.tsdecARroots}, an object of class \code{tsdecARroots} 
returned by \code{type="roots.allocation"}.}
\item{width}{numeric of length two, width of the interval of frequencies allocated to the trend 
and the seasonal components (measured in radians). If a numeric of length one is passed as argument, the same 
width is used for both components.}
%the interval is centered at the fundamental and harmonic seasonal frequencies
\item{min.modulus}{numeric, minimum modulus of the roots assigned to the trend component.}
%\item{n1}{numeric, number of points to be used to draw the contour of the unit circle (by default no contour).}
%\item{n2}{number, number of points to be used to obtain the points of each area.}
\item{xlim}{optional numerics, lower and upper limits of the x-axis.}
\item{ylim}{optional numerics, lower and upper limits of the y-axis.}
\item{units}{character, the units in which the argument of the roots are printed. \code{units="pi"} prints the 
argument in radians as multiples of \eqn{\pi}{pi}.}
\item{digits}{numeric, the number of significant digits to be used by \code{\link[base]{print}}.}
\item{\dots}{further arguments to be passed to \code{\link[graphics]{plot}} or \code{\link[base]{print}}.}
}

\details{
The roots related to cycles with frequency within the range [\eqn{0}, \code{width[1]}] 
are allocated to the trend or transitory component. In particular,if the modulus 
is below \code{min.modulus}, then they are are allocated to the transitory component,
otherwise to the trend.

The seasonal frequencies are defined as \eqn{\omega_j=2\pi j/S}{\omega[j]=2\pi j/S}, 
for \eqn{j=1,...,S-1}, 
where \eqn{S} is the periodicity of the data (e.g., \eqn{S=4} in quarterly data 
and \eqn{S=12} in monthly data).
Roots related to cycles of frequency within the range 
[\eqn{\omega_j - }{\omega[j] - }\code{seasonal.width}, 
\eqn{\omega_j + }{\omega[j] + }\code{seasonal.width}] are
assigned to the seasonal component.
}

\value{
\code{roots.allocation} returns a list of class \code{tsdecARroots}.
\code{plot.tsdecARroots} displays the roots in the complex plane and 
\code{print.tsdecARroots} shows a summary.
}

%\references{}

%\seealso{}

%see give this example, two Real roots (w=0) but one goes to trend and other to transitory component
%set.seed(125)
%y <- arima.sim(n=200, model=list(order=c(2,1,2), ar=c(-0.4,0.2), ma=c(-0.3,0.3)))[-1]
%y <- round(ts(y, frequency=4), 2)
%fit <- arima(y, order=c(2,1,2), fixed=c(-0.33932, 0.26669, -0.45527, 0.32852))

%\examples{}

\keyword{ts}
