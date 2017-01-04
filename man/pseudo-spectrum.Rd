\name{pseudo.spectrum}
\alias{pseudo.spectrum}
\alias{print.tsdecPSP}

\title{Pseudo-Spectrum of an ARIMA Model}

\description{
Compute the polynomials in the numerators of a partial fraction 
decomposition of the pseudo-spectrum in an ARIMA model. 
The polynomials are in terms of the variable \eqn{2\cos\omega}{2*cos(w)}, 
with \eqn{\omega\in [0, 2\pi]}{0 <= w <= 2*pi}.
}

\usage{
pseudo.spectrum(mod, ar)
\method{print}{tsdecPSP}(x, ...)
}

\arguments{
\item{mod}{an object of class \code{Arima}, the fitted model.}
\item{ar}{an object of class \code{tsdecARroots} returned by \code{\link{roots.allocation}}.}
\item{x}{an object of class \code{tsdecPSP} returned by \code{pseudo.spectrum}.}
\item{\dots}{further arguments to be passed to \code{\link[base]{print}}.}
}

\details{
The coefficients of the ARIMA models for each component (e.g., trend, seasonal) 
are obtained from the following relationship.

\deqn{
\sigma^2\frac{\theta(B)\theta(F)}{\phi(B)\phi(F)} = 
\sigma^2_a\frac{\theta_T(B)\theta_T(F)}{\phi_T(B)\phi_T(F)} + 
\sigma^2_b\frac{\theta_S(B)\theta_S(F)}{\phi_S(B)\phi_S(F)} + \sigma^2_e \,,
}{%
sigma^2(theta(B)theta(F))/(phi(B)phi(F)) = 
sigma^2_a(theta_T(B)theta_T(F))/(phi_T(B)phi_T(F)) + 
sigma^2_b(theta_S(B)theta_S(F))(phi_S(B)phi_S(F)) + sigma^2_e,
}

where \eqn{B}{B} is the backshift operator and \eqn{F=B^{-1}}{F=B^(-1)} is the 
forward operator. Each term in the right-hand-side is related to the ARIMA models 
of each one of the unobserved components.

\code{pseudo.spectrum} computes the symmetric polynomials of the type 
\eqn{\varphi(B)\varphi(F)}{varphi(B)varphi(F)} for the polynomials in 
the left-hand-side LHS (based on the fitted model) and for the polynomials 
in the denominators of the right-hand-side RHS 
(based on the allocation of roots of the fitted AR polynomial, 
\code{\link{roots.allocation}}).
Then coefficients in the numerators of the RHS are obtained 
by means of \code{\link{partial.fraction}} .To do so 
the terms in the RHS are multiplied by the denominator in the LHS; 
then, the coefficients of the numerators in the RHS are obtained 
by equating the coefficients of the same order on both sides of 
the relationship (the orders of the unknown polynomials are set 
to one degree lower than those polynomials of the corresponding denominator).
}

\value{
A list of class \code{tsdecPSP} containing:
the quotient of the polynomial division (if the degree of the 
numerator in the LHS is equal or higher than the degree of the denominator);
the coefficients of total polynomials (numerator and denominator in the LHS)
and the denominators in the RHS.
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
\code{\link[stats]{arima}}, \code{\link{partial.fraction}}, 
\code{\link{roots.allocation}}.
}

%\examples{}

\keyword{ts}
