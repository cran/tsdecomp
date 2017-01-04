
dsfilter <- function(x, w, mod, extend = 16)
{
  #this function is the same as f() defined within filtering() below
  #dsfilter() can be used to apply a given filter to other series

  if (extend > 0)
  {
    #if (missing(mod))
    #  stop("missing argument ", sQuote("mod"), 
    #  ", required to extend the series with forecasts")
    startx0 <- start(x)
    endx0 <- end(x)
    pred.right <- stats::predict(mod, n.ahead=extend, se.fit=FALSE)
    revx <- ts(rev(x))
    tsp(revx) <- tsp(x)
    args <- as.list(mod$call)
    args[[1]] <- NULL
    args$x <- revx
    args$fixed <- coef(mod)
    fit2 <- do.call("arima", args)
    pred.left <- stats::predict(fit2, n.ahead=extend, se.fit=FALSE)
    x <- ts(c(pred.left, x, pred.right), 
      frequency=frequency(x), start=time(x)[1]-extend/frequency(x))
  }

  n <- length(x)

  # double-side symmetric filter (from -inf to inf)
  #NOTE do not convolve, just do symmetric with the values from ARMAacov()
  #w <- convolve(w, w, type="open")
  #w <- c(rev(w[-1]), w)
  comp <- ts(stats::filter(c(rep(0, n-1), x, rep(0, n-1)), filter=c(rev(w[-1]), w),
    method="conv", sides=1)[-seq_len(2*(n-1))])
  tsp(comp) <- tsp(x)
  if (extend > 0) {
    window(comp, start=startx0, end=endx0)
  } else 
    comp
}

filtering <- function(x, mod, 
  trend = list(ar=1, ma=1, sigma2=NULL), 
  transitory = list(ar=1, ma=1, sigma2=NULL), 
  seasonal = list(ar=1, ma=1, sigma2=NULL),
  irregular.sigma2 = NULL,
  extend = 16, drift = FALSE)
{
  f <- function(w)
  {
    # same as dsfilter(), "x2" from the scope of filtering() is 
    # used and does not need to be created again

    if (is.null(w))
      return(NULL)
    # double-side symmetric filter (from -inf to inf)

    #NOTE do not convolve, just do symmetric with the values from ARMAacov()
    #w <- convolve(w, w, type="open")
    #w <- c(rev(w[-1]), w)
    comp <- ts(stats::filter(x2, filter=c(rev(w[-1]), w),
      method="conv", sides=1)[-seq_len(nm1x2)])
    tsp(comp) <- tsp(x)
    if (extend > 0) {
      window(comp, start=start(x0), end=end(x0))
    } else 
      comp
  }

  is.trend.linear <- FALSE
  theta <- -mod$model$theta[seq_len(mod$arma[2]+mod$arma[4]*mod$arma[5])]
  # copy of 'mod$model$Delta', it is modified if drift=TRUE
  polyDelta <- c(1, -mod$model$Delta)

  # ARMA coefficients (not including external regressors)
  ncoefs <- sum(mod$arma[seq_len(4)])

  x0 <- x
  if (extend > 0)
  {
    # extend series with forecasts on both ends

    #'drift' is the name given to the external regressor defined in argument 'xreg'
    # include.mean=TRUE is ignoed if D>0, 
    #in that case the drift must be defined through 'xreg'
    
    if ((any(c("intercept", "drift") %in% names(coef(mod)))) && drift)
    {
##FIXME see for models with mod$arma[6]>0 or seasonal differencing

      is.trend.linear <- TRUE
      if ("intercept" %in% names(coef(mod))) {
        linear.trend <- coef(mod)["intercept"]*seq_along(x)
      } else
      if ("drift" %in% names(coef(mod)))
        linear.trend <- coef(mod)["drift"]*seq_along(x)

      x <- x - linear.trend

      mod <- arima(x, order=c(mod$arma[1], mod$arma[6]+1, mod$arma[2]), 
        seasonal=list(order=mod$arma[c(3,7,4)]),
        include.mean=FALSE, fixed=coef(mod)[seq_len(ncoefs)])

      pred.right <- stats::predict(mod, n.ahead=extend, se.fit=FALSE)

    } else  {
##FIXME TODO allow use predictions for external regressors
      newxreg=matrix(0, nrow=extend, ncol=length(coef(mod))-ncoefs)
      pred.right <- stats::predict(mod, n.ahead=extend, se.fit=FALSE, newxreg=newxreg)
    }

    revx <- ts(rev(x))
    tsp(revx) <- tsp(x)
    mod <- arima(revx, order=mod$arma[c(1,6,2)], 
      seasonal=list(order=mod$arma[c(3,7,4)]), 
      include.mean="intercept" %in% names(coef(mod)), 
      fixed=coef(mod)[seq_len(ncoefs)])
    pred.left <- stats::predict(mod, n.ahead=extend, se.fit=FALSE)

    x <- ts(c(rev(pred.left), x, pred.right), 
      frequency=frequency(x), start=time(x)[1]-extend/frequency(x))

  } # end extend > 0

  n <- length(x)
  nm1x2 <- 2*(n-1)
  x2 <- c(rep(0, n-1), x, rep(0, n-1))

  if (!is.null(trend$sigma2)) {
    wtrend <- ARMAacov(ar=theta, ma=polyprod(polyprod(seasonal$ar, transitory$ar), trend$ma)[-1], 
      lag.max=n-1, sigma2=trend$sigma2)
  } else wtrend <- NULL

  if (!is.null(transitory$sigma2)) {
    wtrans <- ARMAacov(ar=theta, ma=polyprod(polyprod(trend$ar, seasonal$ar), transitory$ma)[-1], 
      lag.max=n-1, sigma2=transitory$sigma2)
  } else wtrans <- NULL

  if (!is.null(seasonal$sigma2)) {
    wseas <- ARMAacov(ar=theta, ma=polyprod(polyprod(trend$ar, transitory$ar), seasonal$ma)[-1], 
      lag.max=n-1, sigma2=seasonal$sigma2)
    # filter for the seasonally adjusted series
    wsadj <- c(1-wseas[1], -wseas[-1])
  } else wseas <- wsadj <- NULL

  if (!is.null(irregular.sigma2)) {
    wirreg <- ARMAacov(ar=theta, ma=polyprod(polyDelta, c(1, -mod$model$phi))[-1], 
      lag.max=n-1, sigma2=irregular.sigma2)
  } else wirreg <- NULL

  if (is.trend.linear) {
    if (is.null(wtrend)) {
      comp.trend <- cbind(linear=linear.trend)
    } else
      comp.trend <- cbind(linear=linear.trend, trend=f(wtrend))
  } else
    comp.trend <- f(wtrend)

  res <- structure(
    # series extended on both ends with forecasts
    list(xextended = if(extend > 0) x else NULL,
    # one-side filters
    filters = cbind(trend=wtrend, transitory=wtrans, seasonal=wseas, 
      sadj=wsadj, irregular=wirreg),
    # filtered components
    components=cbind(observed=x0, 
      trend=comp.trend, transitory=f(wtrans), seasonal=f(wseas), 
      sadj=f(wsadj), irregular=f(wirreg))), 
    class="tsdecFilter")

##NOTE this is required because 
#concatenation of 'ts' with NULL does not keep the correct column names
#(this has no effect if there is no NULL elements)  
  if (is.trend.linear) 
  {
    if (NCOL(comp.trend) > 1) {
      colnames(res$components) <- c("observed", "linear.trend", colnames(res$filters))
    } else
      colnames(res$components) <- c("observed", "trend", colnames(res$filters))
  } else
    colnames(res$components) <- c("observed", colnames(res$filters))

  res
}
