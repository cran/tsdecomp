\name{canonical.decomposition}
\alias{canonical.decomposition}
\alias{print.tsdecCanDec}

\title{Canonical Decomposition}

\description{
Given the partial fraction decomposition of the pseudo-spectrum, 
the canonical decomposition allocates the variance of each component 
so that the variance of the irregular is maximised. 
Then, the coefficients of the numerators in the pseudo-spectrum 
(relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}) 
are converted into the MA coefficients of the model for each component 
by means of \code{\link[tsdecomp]{acgf2poly}}.
}

\usage{
canonical.decomposition(num.trend, den.trend, 
  num.trans, den.trans, num.seas, den.seas, quotient, optim.tol = 1e-04, ...)
\method{print}{tsdecCanDec}(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
}

\arguments{
\item{num.trend}{numeric vector, the coefficients of the MA polynomial related to 
the trend component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{den.trend}{numeric vector, the coefficients of the AR polynomial related to 
the trend component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{num.trans}{numeric vector, the coefficients of the MA polynomial related to 
the transitory component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{den.trans}{numeric vector, the coefficients of the AR polynomial related to 
the transitory component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{num.seas}{numeric vector, the coefficients of the MA polynomial related to 
the seasonal component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{den.seas}{numeric vector, the coefficients of the AR polynomial related to 
the seasonal component in the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.}
\item{quotient}{numeric vector, the quotient of the polynomial division of the polynomials 
in the LHS of the relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}.
(Different from zero only when 
the degree of the MA polynomial is equal or higher than the degree of the AR polynomial 
in the fitted model).}
\item{optim.tol}{numeric, the convergence tolerance to be used by \code{\link[stats]{optimize}}.}
\item{units}{character, the units in which the argument of the roots are printed. \code{units="pi"} prints the 
argument in radians as multiples of \eqn{\pi}{pi}.}
\item{x}{an object of class \code{tsdecCanDec} returned by \code{canonical.decomposition}.}
\item{digits}{numeric, the number of significant digits to be used by \code{\link[base]{print}}.}
\item{\dots}{further arguments to be passed to \code{\link[tsdecomp]{poly2acgf}} or 
\code{\link[base]{print}}.}
}

%\details{}

\value{
An object of class \code{tsdecCanDec} containing 
the MA coefficients of the ARIMA models obtained for the unobserved components
(e.g., trend, seasonal) and the variance of the corresponding disturbance terms.
}

\references{
Burman, J. P. (1980) 
\sQuote{Seasonal Adjustment by Signal Extraction}.
\emph{Journal of the Royal Statistical Society. Series A (General)},
\bold{143}(3), pp. 321-337.
\doi{10.2307/2982132}.

Hillmer, S. C. and Tiao, G. C. (1982) 
\sQuote{An ARIMA-Model-Based Approach to Seasonal Adjustment}.
\emph{Journal of the American Statistical Association},
\bold{77}(377), pp. 63-70.
\doi{10.1080/01621459.1982.10477767}.
}

\seealso{
\code{\link[tsdecomp]{acgf2poly}}, \code{\link[tsdecomp]{pseudo.spectrum}}, 
\code{\link[stats]{optimize}}.
}

%\examples{}

\keyword{ts}
