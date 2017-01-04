\name{partial.fraction}
\alias{partial.fraction}

\title{Partial Fraction Decomposition}

\description{
Partial fraction decomposition of the pseudo-spectrum of a fitted ARIMA model.
}

\usage{
partial.fraction(numerator, den.trend, den.transitory, den.seasonal)
}

\arguments{
\item{numerator}{numeric vector containing the coefficients of the numerator of the 
ratio of polynomials to be decomposed (numerator in the left-hand-side of the 
relationship given in \code{\link[tsdecomp]{pseudo.spectrum}}).}
\item{den.trend}{numeric vector containing the coefficients of the denominator in 
the partial fraction related to the trend component.}
\item{den.transitory}{numeric vector containing the coefficients of the denominator in 
the partial fraction related to the transitory component.}
\item{den.seasonal}{numeric vector containing the coefficients of the denominator in 
the partial fraction related to the seasonal component.}
}

%\details{
%The relationship in which the decomposition relies is given in 
%\code{\link[tsdecomp]{pseudo.spectrum}}.
%}

\value{
A list containing the system of equations which is solved and 
the numerators of the partial fractions related, respectively, 
to the trend, transitory and seasonal components.
}

%\references{}

\seealso{
\code{\link[tsdecomp]{pseudo.spectrum}}.
}

%\examples{}

\keyword{ts}
