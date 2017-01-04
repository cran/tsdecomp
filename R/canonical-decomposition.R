
canonical.decomposition <- function(num.trend, den.trend, 
  num.trans, den.trans, num.seas, den.seas, quotient, optim.tol = 1e-04, ...)
{
  f <- function(x, p1, p2) polyeval(p1, x) / polyeval(p2, x)

  trend <- trans <- seas <- list(coef=1) #list(roots=NULL, coef=NULL, sigma2=NULL)
  trend.minval <- trans.minval <- seas.minval <- 0

  if (!is.null(num.trend))
  {
    res <- stats::optimize(f=f, interval=c(-2, 2), p1=num.trend, p2=den.trend, tol=optim.tol)
    trend.minval <- min(res$obj, f(-2, num.trend, den.trend))
    num.trend <- c(num.trend, 0) - den.trend*trend.minval
    trend <- poly2acgf(num.trend, ...)
  }

  if (!is.null(num.trans))
  {
    trans.minval <- stats::optimize(f=f, interval=c(-2, 2), 
      p1=num.trans, p2=den.trans, tol=optim.tol)$obj
    num.trans <- c(num.trans, 0) - den.trans*trans.minval
    trans <- poly2acgf(num.trans, ...)
  }

  if (!is.null(num.seas))
  {
    res <- stats::optimize(f=f, interval=c(-2, 2), p1=num.seas, p2=den.seas, tol=optim.tol)
    seas.minval <- min(res$obj, f(2, num.seas, den.seas))
    num.seas <- c(num.seas, 0) - den.seas*seas.minval
    seas <- poly2acgf(num.seas, ...)
  }

  if ((lq <- length(quotient)) > 1)
  {
    # do the canonical decomposition of the term below and 
    # assign 'irregular.minval' to the variance of the irregular component
    irregular <- quotient + 
      c(trend.minval + trans.minval + seas.minval, rep(0, lq-1))

    irregular.minval <- min(
      stats::optimize(f=polyeval, interval=c(-2, 2), p=irregular, tol=optim.tol)$obj, 
      polyeval(p=irregular, x=-2), polyeval(p=irregular, x=2))

    irregular <- irregular - c(irregular.minval, rep(0, lq-1))
    tmp <- poly2acgf(irregular, ...)

    if (is.null(num.trans))
    {
      # assign to a transitory component
      cd.num.trans <- tmp$coef
      sigma2.trans <- tmp$sigma2
    } else {
      # I would have expected that the AR roots allocated to the transitory 
      # would have captured this remaining component as well;
      # If this happens, see create explicity a model for the irregular or 
      # add up both transitory components
      #cd.num.irregular <- tmp$coef
      #irregular.sigma2 <- tmp$sigma2
      #stop("a transitory component showed up both in AR allocation and partial fraction")
      stop("unsupported model: two transitory components were found")
    }

    irregular.sigma2 <- irregular.minval

  } else # length(quotient) == 1 (degree of numerator <= degree of denominator)
    irregular.sigma2 <- quotient + trend.minval + trans.minval + seas.minval

  # recover coefficients from spectrum to ACGF of the MA part 
  # (it belongs to the MA part because the AR part vanished when multiplying 
  # by 'den.trend' and 'den.seas' when 'num.trend' and 'num.seas' updated above)

  if (irregular.sigma2 < 0)
  {
##FIXME in theory this may happen for some models, 
#but check the procedure if this happens
    stop("decomposition not admissible (negative spectrum of the irregular component)")
  }

  structure(
    list(trend=trend, transitory=trans, seasonal=seas, irregular.sigma2=irregular.sigma2),
    class="tsdecCanDec")
}

print.tsdecCanDec <- function(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
{
  units <- match.arg(units)

  # MA polynomials

  cat(paste("MA polynomials\n--------------\n\n"))

  if (!is.null(x$trend$sigma2))
  {
    cat(paste0("  Trend:\n  ", 
      polystring(x$trend$coef, varchar="L", brackets=TRUE, ndec=3), 
      "a_t,  a_t ~ IID(0, ", round(x$trend$sigma2, digits=digits), ")\n"))
  }

  if (!is.null(x$transitory$sigma2))
  {
    cat(paste0("  Transitory:\n  ", 
      polystring(x$transitory$coef, varchar="L", brackets=TRUE, ndec=3), 
      "b_t,  b_t ~ IID(0, ", round(x$transitory$sigma2, digits=digits), ")\n"))
  }

  if (!is.null(x$seasonal$sigma2))
  {
    cat(paste0("  Seasonal:\n  ", 
      polystring(x$seasonal$coef, varchar="L", brackets=TRUE, ndec=3), 
      "c_t,  c_t ~ IID(0, ", round(x$seasonal$sigma2, digits=digits), ")\n"))
  }

  # Variances

  cat(paste("\nVariances\n---------\n\n"))
  print(c(trend=x$trend$sigma2, transitory=x$transitory$sigma2, 
    seasonal=x$seasonal$sigma2, irregular=x$irregular.sigma2), digits=digits, ...)

  # Roots

  cat(paste("\nRoots\n-----\n\n"))
  m1 <- if (inherits(x$trend, "tsdecMAroots"))
    print(x$trend, units=units, echo=FALSE) else NULL
  m2 <- if (inherits(x$transitory, "tsdecMAroots"))
    print(x$transitory, units=units, echo=FALSE, ...) else NULL
  m3 <- if (inherits(x$seasonal, "tsdecMAroots"))
    print(x$seasonal, units=units, echo=FALSE, ...) else NULL
  m <- cbind(Component=c(rep("trend", max(0, nrow(m1))), 
    rep("transitory", max(0, nrow(m2))), rep("seasonal", max(0, nrow(m3)))), 
    rbind(m1,m2,m3))

  print(m, digits=digits, ...) #right=FALSE
}
