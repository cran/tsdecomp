\name{tsdecomp-package}
\docType{package}
\alias{tsdecomp-package}
\alias{tsdecomp}

\title{ARIMA-Model-Based Decomposition of Time Series Data}

\description{
ARIMA-model-based decomposition of a time series.
}

\details{
The methods implemented in the package are developed and described (among others)
in the references given below. 
The package is mainly intended for annual, quarterly and monthly time series.
The bottom line of the procedure can be summarized as follows.
An ARIMA model is fitted to the observed series. Then the pseudo-spectrum 
of the model is computed and decomposed into partial fractions. 
This gives ARIMA models for the unobserved components 
(e.g., trend, seasonal and irregular), 
which are then used to obtain the weights of double-sided linear filters 
upon which estimates of the components are obtained.

For practical purposes, the main function provided in the package is 
\code{\link[tsdecomp]{ARIMAdec}}. This function relies on other 
procedures that implement different stages of the procedure:
\code{\link[tsdecomp]{roots.allocation}},
\code{\link[tsdecomp]{pseudo.spectrum}},
\code{\link[tsdecomp]{partial.fraction}},
\code{\link[tsdecomp]{canonical.decomposition}},
\code{\link[tsdecomp]{filtering}}.

An introduction to the methodology and the package 
in the form of a vignette is available here:

\url{https://www.jalobe.com/doc/tsdecomp.pdf}
}

\author{
Javier López-de-Lacalle \email{javlacalle@yahoo.es}

\url{https://jalobe.com}
}

\references{
Box, G. E. P., Hillmer, S. C. and Tiao, G. C. (1978)
\sQuote{Analysis and Modeling of Seasonal Time Series} in 
\emph{Seasonal Analysis of Economic Time Series}, Editor Zellner, A.
pp. 309-334. U.S. Dept. of Commerce - Bureau of the Census.
\url{http://www.nber.org/chapters/c3904.pdf}

Brockwell, P. J. and Davis, R. A. (1991) 
\emph{Time Series: Theory and Methods}, 
Second Edition. Springer.
\doi{10.1007/978-1-4419-0320-4}

Burman, J. P. (1980) 
\sQuote{Seasonal Adjustment by Signal Extraction}.
\emph{Journal of the Royal Statistical Society. Series A (General)},
\bold{143}(3), pp. 321-337.
\doi{10.2307/2982132}

Gómez, V. and Maravall, A. (2001)
\sQuote{Programs \pkg{TRAMO} and \pkg{SEATS}. Instructions for the User (Beta Version: June 1997)}.
Ministerio de Economía y Hacienda. Dirección General de Análisis y Programación Presupuestaria,
Working paper SGAPE-97001.
\url{http://www.bde.es/f/webbde/SES/servicio/Programas_estadisticos_y_econometricos/Programas/ficheros/manualdos.pdf}

Gómez, V. (2015) 
\sQuote{\pkg{SSMMATLAB}: A Set of MATLAB Programs for the Statistical Analysis of State Space Models}.
\emph{Journal of Statistical Software}, \bold{66}(1), pp. 1-37.
\doi{10.18637/jss.v066.i09}
\url{http://www.sepg.pap.minhap.gob.es/sitios/sepg/en-GB/Presupuestos/Documentacion/Paginas/SSMMATLAB.aspx}

Hillmer, S. C. and Tiao, G. C. (1982) 
\sQuote{An ARIMA-Model-Based Approach to Seasonal Adjustment}.
\emph{Journal of the American Statistical Association},
\bold{77}(377), pp. 63-70.
\doi{10.1080/01621459.1982.10477767}

Maravall, A. and Pierce, D. A. (1987)
\sQuote{A Prototypical Seasonal Adjustment Model}.
\emph{Journal of Time Series Analysis}, \bold{8}(2), pp.177-193.
\doi{10.1111/j.1467-9892.1987.tb00431.x}

Planas, C. (1997) 
\emph{Applied Time Series Analysis: Modelling, Forecasting, Unobserved Components Analysis and the Wiener-Kolmogorov Filter}.
Eurostat: Series E, Methods.
\url{https://bookshop.europa.eu/en/applied-time-series-analysis-pbCA0897484/}

Pollock, D. S. G. (1999)
\emph{A Handbook of Time-Series Analysis Signal Processing and Dynamics}. 
Academic Press. 
\doi{10.1016/B978-012560990-6/50002-6}
}

\keyword{package, ts}
