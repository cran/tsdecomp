\name{acgf2poly}
\alias{acgf2poly}
\alias{poly2acgf}
\alias{print.tsdecMAroots}

\title{Change of Variable in the AutoCovariance Generating Function}

\description{
Change of variable in the autocovariance generating function (ACGF). 
This transformation allows the pseudo-spectrum to be represented as 
a polynomial liable to be decomposed in partial fractions.
}

\usage{
acgf2poly(x)
poly2acgf(x, type=c("roots2poly", "acov2ma"), tol = 1e-16, maxiter = 100, 
  init.tol=1e-05, init.maxiter=100)
\method{print}{tsdecMAroots}(x, units = c("radians", "degrees", "pi"), digits = 4, echo = TRUE, ...)
}

\arguments{
\item{x}{numeric vector of autocovariances;
for \code{poly2acgf}, an object of class \code{tsdecMAroots} returned by \code{type="roots2poly"}}
\item{type}{character string selecting the method to undo the transformation.}
\item{tol}{convergence tolerance to be used by \code{\link{acov2ma}}.}
\item{maxiter}{maximum number of iterations allowed in \code{\link{acov2ma}}.}
\item{init.tol}{convergence tolerance to be used by \code{\link{acov2ma.init}}.}
\item{init.maxiter}{maximum number of iterations allowed in \code{\link{acov2ma.init}}.}
\item{units}{character, the units in which the argument of the roots are printed. \code{units="pi"} prints the 
argument in radians as multiples of \eqn{\pi}{pi}.}
\item{digits}{numeric, the number of significant digits to be used by \code{\link[base]{print}}.}
\item{echo}{logical, if \code{TRUE} the output is printed, otherwise a \code{\link[base]{invisible}} 
copy of the matrix summarizing the results obtained by \code{poly2acgf} is returned.}
\item{\dots}{further arguments to be passed to \code{\link[base]{print}}.}
}

\details{
The ACGF is defined as a power series where the coefficients are the autocovariances 
\eqn{\gamma_\tau}{g[i]}:
%
%\deqn{
%\gamma(L) = \sum_{\tau=-\infty}^\infty \gamma_\tau L^\tau}{%
%g(L) = \sum_{i=-infty}^infty g[i] L^i}
%
%It can also be expressed as:
%
\deqn{
\gamma(z) = \gamma_0 + \gamma_1(z+z^{-1}) + \gamma_2(z^2+z^{-2}) + \gamma_3(z^3+z^{-3}) + \cdots}{%
g(z) = g[0] + g[1](z+z^(-1)) + g[2](z^2+z^(-2)) + g[3](z^3+z^(-3)) + ...}
%
where \eqn{z}{z} is a complex variable.

Replacing \eqn{z}{z} by \eqn{e^{-i\omega}}{e^(-i*w)} with \eqn{\omega\in[0,2\pi]}{0 <= w <= 2*pi}
yields the spectral density multiplied by \eqn{2\pi}{2*pi}.
This gives a power series in the variable \eqn{2\cos(\omega j)}{2*cos(w*j)}
(note that for \eqn{z=e^{-i\omega}}{z=e^(-i*w)}, which has unit modulus, 
the inverse \eqn{1/z}{1/z} is the complex-conjugate of \eqn{z}{z}):
%
\deqn{z^j + z^{-j} = \cos(\omega j) + i\sin(\omega j) + \cos(\omega j) - i\sin(\omega j) = 
2\cos(\omega j)\,.}{%
z^j + z^(-j) = cos(w*j) + i*sin(w*j) + cos(w*j) - i*sin(w*j) = 2*cos(w*j).}

\code{acgf2poly} transforms the following expression in the variable \eqn{2\cos(\omega j)}{2*cos(w*j)}:
%
\deqn{A(2\cos(j\omega)) = a_0 + a_1 2\cos(\omega) + a_2 2\cos(2\omega) + \cdots + a_n 2\cos(n\omega)}{%
A(2*cos(j*w)) = a[0] + a[1] * 2*cos(w) + a[2] 2*cos(2*w) + ... + a[n] 2*cos(n*w)}
%
into a polynomial in the variable \eqn{x=2\cos(\omega)}{x=2*cos(w)}:
%
\deqn{
B(x) = b_0 + b_1 x + b_2 x^2 + \cdots + b_n x^n\,.}{%
B(x) = b[0] + b[1]*x + b[2]*x^2 + ... + b[n] x^n.}

\code{poly2acgf} recovers the vector of autocovariances by 
undoing the above transformation and computes the coefficients 
and the variance of the innovations of the moving average model 
related to those autocovariances. 
Two methods can be employed.
1) \code{type="acov2ma"}: this method recovers the autocovariances 
by undoing the change of variable; then, the the autocovariances 
are converted to the coefficients of a moving average 
by means of \code{\link{acov2ma}}. 
In the presence of non-invertible roots, 
this method may experience difficulties to converge.
%
2) \code{type="roots2poly"}: this method does not explicitly undo 
the change of variable \code{acgf2poly} (i.e., the vector of 
autocovariances is not recovered). Instead, 
the roots of the moving average polynomial \eqn{\theta(L)}{theta(L)}
are obtained from the polynomial \eqn{\theta(L)\theta(L^{-1})}{theta(L)*theta(L^(-1))}, 
where the coefficients are in terms of the polynomial \eqn{B(x)}{B(x)} 
defined above; then, the coefficients of the moving average model are 
computed by means of \code{\link{roots2poly}}.
}

\value{
\code{acgf2poly} returns the transformed vector of coefficients.

\code{poly2acgf} returns an object of class \code{tsdecMAroots}
containing the coefficients and the variance of the 
innovations in the moving average model related to the autocovariances 
underlying the transformed coefficients.
\code{print.tsdecMAroots} prints a summary of the results computed 
by \code{poly2acgf}.
}

\note{
Method \code{type="roots2poly"} in \code{poly2acgf} is based on 
algorithm \code{pu2ma} in the software SSMMATLAB by Gómez, V. 
URL 
\url{http://www.sepg.pap.minhap.gob.es/sitios/sepg/en-GB/Presupuestos/Documentacion/Paginas/SSMMATLAB.aspx}.
}

%\references{}

\seealso{
\code{\link{acov2ma}}, \code{\link{roots2poly}}.
}

\examples{
# the matrix 'm' performs the mapping from the original 
# to the transformed coefficients
n <- 30
m <- diag(1, n, n)
n2 <- n - 2
j <- -1
tmp <- seq.int(2, n-1)
for (i in seq.int(3, n-2, 2))
{
  id <- cbind(seq_len(n2),seq.int(i,n))
  m[id] <- j * tmp
  n2 <- n2 - 2
  j <- -1 * j
  tmp <- cumsum(tmp[seq_len(n2)])
}
if (2*floor(n/2) == n) {  # if (n %% 2 == 0)
  m[cbind(seq_len(n2),seq.int(n-1,n))] <- j * tmp 
} else 
  m[1,n] <- j * tmp
m[1:10,1:10]

# equivalence of the original and transformed coefficients,
# example with an ARMA(2,1) model
#
# method 1: compute the spectral density upon the 
# the theoretical autocovariances ('gamma') of the ARMA model
gamma <- ARMAacov(ar=c(0.8,-0.6), ma=0.4, lag.max=n-1)
w <- seq(0, pi, len=length(gamma))
spec1 <- rep(gamma[1], length(w))
for (i in seq_along(w))
{
  z <- 2*cos(w[i] * seq_len(length(gamma)-1))
  spec1[i] <- spec1[i] + sum(gamma[seq.int(2, n)] * z)
}
spec1 <- spec1/(2*pi)
#plot(w, spec1)

# method 2: compute the spectral density upon the 
# transformed coefficients
newcoefs <- m %*% gamma
spec2 <- rep(newcoefs[1], length(w))
for (i in seq_along(w))
{
  x <- (2*cos(w[i]))^seq_len(n-1)
  spec2[i] <- spec2[i] + sum(newcoefs[seq.int(2, n)] * x)
}
spec2 <- spec2/(2*pi)

# both representations are equivalent
all.equal(spec1, spec2, check.names=FALSE)
#[1] TRUE

# the original coefficients (the autocovariances) 
# can be recovered premultiplying by the inverse of the 
# transformation matrix 'm'
all.equal(c(solve(m) \%*\% newcoefs), gamma, check.names=FALSE)
#[1] TRUE
}

\keyword{ts}
