
##FIXME TODO irregular component
##FIXME check with more models (see examples in documentation)

compare.acf <- function(x, mod, lag.max = 12, ...)
{
  do.trend <- !is.null(x$ma$trend$sigma2)
  do.trans <- !is.null(x$ma$transitory$sigma2)
  do.seas <- !is.null(x$ma$seasonal$sigma2)

  # ACF of the ARMA models for the components
  # theoretical ACF, assuming 'mod' is the 'true' or correct
  # model for the data

  a1 <- a2 <- a3 <- NULL
  if (do.trend)
    a1 <- stats::ARMAacf(ar=-x$ar$polys.stationary$trend[-1], 
      ma=x$ma$trend$coef[-1], lag.max=lag.max)[-1]
  if (do.trans)
    a2 <- stats::ARMAacf(ar=-x$ar$polys.stationary$transitory[-1], 
      ma=x$ma$transitory$coef[-1], lag.max=lag.max)[-1]
  if (do.seas)
    a3 <- stats::ARMAacf(ar=-x$ar$polys.stationary$seasonal[-1], 
      ma=x$ma$seasonal$coef[-1], lag.max=lag.max)[-1]

  # ACF of the estimators for the components

  b1 <- b2 <- b3 <- NULL
  theta <- c(1, mod$model$theta[seq_len(mod$arma[2]+mod$arma[4]*mod$arma[5])])
  #theta <- roots2poly(1/polyroot(theta))
  r <- polyroot(theta)
  id <- which(Mod(r) >= 1)
  theta <- if (length(id) > 0) -roots2poly(r[-id])[-1]  else 
    roots2poly(1/polyroot(theta))
  if (length(theta) == 0)
    theta <- 1

  if (do.trend)
  {
    ar.coefs <- -polyprod(x$ar$polys.stationary$trend, theta)[-1]
    ma.coefs <- stats::convolve(x$ma$trend$coef, x$ma$trend$coef, type="open") 
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$transitory)))
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$seasonal)))
    b1 <- stats::ARMAacf(ar=ar.coefs, ma=ma.coefs[-1], lag.max=lag.max)[-1]
    #b <- ARMAacov(ar=ar.coefs, ma=ma.coefs[-1], lag.max=lag.max, sigma2=x$ma$trend$sigma2)
    #b/b[1]
  }

  if (do.trans)
  {
    ar.coefs <- -polyprod(x$ar$transitory, theta)[-1]
    ma.coefs <- stats::convolve(x$ma$transitory$coef, x$ma$transitory$coef, type="open") 
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$trend)))
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$seasonal)))
    b2 <- stats::ARMAacf(ar=ar.coefs, ma=ma.coefs[-1], lag.max=lag.max)[-1]
  }

  if (do.seas)
  {
    ar.coefs <- -polyprod(x$ar$polys.stationary$seasonal, theta)[-1]
    ma.coefs <- stats::convolve(x$ma$seasonal$coef, x$ma$seasonal$coef, type="open") 
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$trend)))
    ma.coefs <- polyprod(ma.coefs, roots2poly(1/polyroot(x$ar$transitory)))
    b3 <- stats::ARMAacf(ar=ar.coefs, ma=ma.coefs[-1], lag.max=lag.max)[-1]
  }

  # ACF of the estimates, empirical/filtered components

  c1 <- c2 <- c3 <- NULL
  if (do.trend)
  {
    f <- x$ar$polys.nonstat$trend
    if (!is.null(x$ar$polys.nonstat$trend) && !identical(f, 1)) 
    {
      tmp <- stats::filter(x$components[,"trend"], filter=f, method="conv", sides=1)
      tmp <- tmp[-seq_along(f[-1])] # remove NAs before calling to 'acf'
    } else 
      tmp <- x$components[,"trend"]
    c1 <- stats::acf(tmp, lag.max=lag.max, plot=FALSE, ...)$acf[-1,,1]
  }

  if (do.trans)
    c2 <- stats::acf(x$components[,"transitory"], lag.max=lag.max, plot=FALSE, ...)$acf[-1,,1]

  if (do.seas)
  {
    f <- x$ar$polys.nonstat$seasonal  
    if (!is.null(x$ar$polys.nonstat$seasonal) && !identical(f, 1)) 
    {
      tmp <- stats::filter(x$components[,"seasonal"], filter=f, method="conv", sides=1)
      tmp <- tmp[-seq_along(f[-1])] # remove NAs before calling to 'acf'
    } else 
      tmp <- x$components[,"seasonal"]
    c3 <- stats::acf(tmp, lag.max=lag.max, plot=FALSE, ...)$acf[-1,,1]
  }

  c123 <- cbind(trend=c1, transitory=c2, seasonal=c3)
  rownames(c123) <- seq_len(nrow(c123))

  structure(list(nobs=mod$nobs,
    theoretical=cbind(trend=a1, transitory=a2, seasonal=a3),
    estimator=cbind(trend=b1, transitory=b2, seasonal=b3),
    empirical=c123), 
    class="tsdecAcf")
}

plot.tsdecAcf <- function(x, component = c("trend", "transitory", "seasonal"), ci = 0.95, 
  ci.type = c("ma", "white"), ci.class = c("estimator", "theoretical", "empirical"), 
  plot = TRUE, ...)
{
  component <- match.arg(component)
  if (!(component %in% colnames(x$theoretical)))
    stop(paste("the model does not contain component", sQuote(component)))

  ci.type <- match.arg(ci.type)
  ci.class <- match.arg(ci.class)
  lag.max <- nrow(x$theoretical)

  if (ci > 0)
  {
    # confidence bands
    # based on stats::plot.acf

    ci0 <- stats::qnorm((1 + ci)/2) / sqrt(x$nobs)
    #argument 'ci.class' is ignored if ci.class=='white'

    if (ci.type == "ma") # Bartlett approximation
    {
      CI <- switch(ci.class, 
        "estimator" = ci0 * sqrt(cumsum(c(1, 2*x$estimator[-lag.max,component]^2))),
        "theoretical" = ci0 * sqrt(cumsum(c(1, 2*x$theoretical[-lag.max,component]^2))), 
        "empirical" = ci0 * sqrt(cumsum(c(1, 2*x$empirical[-lag.max,component]^2))))
      names(CI) <- seq_len(lag.max)

    } else
      CI <- ci0

    ylim <- range(c(-CI, CI, x$theoretical[,component], 
      x$estimator[,component], x$empirical[,component]))

    if (ci.type == "ma") {
      ylim[1] <- if (sign(ylim[1]) == 1) max(-1, ylim[1]*0.9) else max(-1, ylim[1]*1.1)
      ylim[2] <- if (sign(ylim[2]) == 1) min(1, ylim[2]*1.1) else min(1, ylim[2]*0.9)
    } else {
      ylim[1] <- if (sign(ylim[1]) == 1) max(-1, ylim[1]*0.8) else max(-1, ylim[1]*1.2)
      ylim[2] <- if (sign(ylim[2]) == 1) min(1, ylim[2]*1.2) else min(1, ylim[2]*0.8)
    }

  } else {
    CI <- NULL
    ylim <- NULL
  }

  if (plot)
  {
    plot(x$theoretical[,component], type="p", xlab="Lag", ylab="", xaxt="n", 
      main=paste0("ACF of ", component), ylim=ylim, ...)
    axis(side=1, at=seq_len(lag.max), labels=seq_len(lag.max)) #labels=seq.int(0, lag.max-1)
    lines(x$estimator[,component], type="h", col="blue")
    lines(x$empirical[,component], type="p", col="red")
    graphics::legend("topright", legend=c("theoretical", "estimator", "empirical"), 
      col=c("black", "blue", "red"), lty=c(0,1,0), pch=c(1,-1,1), bty="n")
    tmp <- if (sign(x$theoretical[1,component]) == 1) "bottomleft" else "bottomright"
    graphics::legend(tmp, legend=paste0(ci*100, "% confidence bands (", ci.class, ")"), 
      col=c("red", "red"), lty=c(2,2), bty="n")
    abline(h=0, lty=2)

    # confidence bands

    if (ci > 0)
    {
      if (ci.type == "white")
      {
        abline(h=c(CI, -CI), col="red", lty=2)
      } else { # ci.type == "ma"
        lines(CI, col="red", lty=2)
        lines(-CI, col="red", lty=2)
      }
    }
  }

  if (ci > 0)
  {
    tmp <- switch(ci.class, "estimator" = x$estimator[,component],
      "theoretical" = x$theoretical[,component], "empirical" = x$empirical[,component])
    CI <- cbind(-CI, tmp, CI)
    colnames(CI) <- c(paste0(ci*100, "% lower"), component, paste0(ci*100, "% upper"))
  }

  invisible(CI)
}
