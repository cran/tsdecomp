\name{polyeval}
\alias{polydiv}
\alias{polyeval}
\alias{polyprod}
\alias{polystring}
\alias{roots2poly}

\title{Polynomial Operations and Utilities}

\description{
Polynomial operations and utilities.
}

\usage{
polystring(x, varchar = "x", brackets = FALSE, ndec = 2, emptychar = "")
polyeval(p, x)
polyprod(x, y, tol = 1.490116e-08) 
polydiv(x, y)
roots2poly(x)
}

\arguments{
\item{x}{numeric vector containing the coefficients of the polynomial (in increasing order and 
without gaps). For \code{polyeval}, this is the point at which the polynomial is to be evaluated.
For \code{roots2poly}, this is a numeric vector containing the roots of the polynomial.}
\item{y}{numeric vector containing the coefficients of the polynomial (in increasing order and 
without gaps).}
\item{p}{numeric vector containing the coefficients of the polynomial (in increasing order and 
without gaps).}
\item{varchar}{character string, the label to be printed representing the variable of the polynomial, 
defaults to \code{"x"}.}
\item{brackets}{logical, if \code{TRUE} the polynomial is printed within parentheses.}
\item{ndec}{integer, coefficients are rounded up to this number of decimals.}
\item{emptychar}{the character string to be printed if the polynomial is empty.}
\item{tol}{a numeric, tolerance below which coefficients are set to zero.}
}

\details{
\code{polystring} returns a string of a numeric vector in the format of a polynomial.

\code{polyeval} evaluates the polynomial defined in the vector of coefficients \code{p} 
at the point \code{x}.

\code{polyprod} performs polynomial multiplication.

\code{polydiv} performs polynomial division (returning the quotient and the remainder).

\code{roots2poly} computes the coefficients of a polynomial from its roots.
}

%\value{}

\note{
\code{polyprod} is based on \code{\link[stats]{convolve}}; it is equivalent to
\code{convolve(x, rev(y), type="open")}.

\code{roots2poly} is based on \code{poly.from.zeros()} in package \CRANpkg{polynom}.
}

%\references{}

\seealso{
\code{\link[base]{polyroot}}.
}

\examples{
# print a fitted ARMA model
set.seed(123)
y <- arima.sim(n=120, model=list(ar=c(0.8, -0.3), ma=0.6))
fit <- arima(y, order=c(2,0,1), include.mean=FALSE)
cat(paste0(
polystring(c(1, -fit$model$phi), brackets=TRUE, ndec=3), "y_t = ",
polystring(c(1, fit$model$theta), brackets=TRUE, ndec=3), "e_t\n"))

# convert roots to coefficients
p <- c(1, 0.8, -0.3)
cat(polystring(p))
r <- polyroot(p)
roots2poly(r)
}

\keyword{ts}
