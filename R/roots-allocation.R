
roots.allocation <- function(x, width = c(0.035, 0.035), min.modulus = 0.4)
{
  arcoefs <- -coef(x)[seq_len(x$arma[1])]
  sarcoefs <- -coef(x)[seq.int(x$arma[1]+x$arma[2]+1, length.out=x$arma[3])]
  S <- x$arma[5]
  d <- x$arma[6]
  D <- x$arma[7]

  if (length(width) == 1)
    width <- c(width, width)

  wseas <- seq_len(S-1)*2*pi/S

  # trend-cycle component, non-stationary part 

  if (d+D > 0)
  {
    rtrend.ns <- rep(1, d+D)

    ptrend.ns <- if (d+D > 1)
      roots2poly(rep(1, d+D)) else c(1, -1)

  } else {
    rtrend.ns <- NULL
    ptrend.ns <- 1
  }

  # seasonal component, non-stationary part 

  if (D > 0)
  {
    rseas.ns <- rep(polyroot(rep(1, S)), D)

    tmp <- pseas.ns <- rep(1, S)
    for (i in seq_len(D-1))
      pseas.ns <- stats::convolve(pseas.ns, tmp, type="open")
      #pseas.ns <- polyprod(pseas.ns, tmp)
      #pseas.ns <- .Call(stats:::C_TSconv, pseas.ns, tmp)
  } else  {
    rseas.ns <- NULL
    pseas.ns <- 1
  }

  # AR polynomial, stationary part 

  if (length(arcoefs) > 0)
  {
    r <- polyroot(c(1, arcoefs))

    if (S > 1)
    {
      tmp <- findInterval(abs(Arg(r)), c(rbind(c(0, wseas - width[2]), 
        c(width[1], wseas + width[2]))), rightmost.closed=TRUE)
      idaux <- which(tmp == 1)
      id.tc <- idaux[which(1/Mod(r[idaux]) > min.modulus)]
      id.trans <- idaux[which(1/Mod(r[idaux]) <= min.modulus)]
      tmp[id.tc] <- NA
      tmp[id.trans] <- NA
      tmp <- tmp %% 2
      id.trans <- c(id.trans, which(tmp == 0))
      id.seas <- which(tmp == 1)

#debug
stopifnot(!any(duplicated(c(id.tc, id.trans, id.seas))))

    } else # S == 1
    {
      # if S=1, allocate positive real roots to trend and remaining to transitory
      tmp <- which(Im(r) == 0)
      id.tc <- intersect(tmp, which(Re(r) > 0))
      id.trans <- intersect(tmp, which(Re(r) < 0))
      id.trans <- c(id.trans, which(Im(r) != 0))
      id.seas <- integer(0)
    }

    rtrend.s <- r[id.tc]
    rtrans <- r[id.trans]
    rseas.s <- r[id.seas]

  } else {
    rtrend.s <- rtrans <- rseas.s <- NULL
  }

  # seasonal AR polynomial, stationary part 

  #same code as above (except for 'r')
  if (length(sarcoefs) > 0)
  {
    #alternatively define for example non-seasonal model of order S, e.g. ar=c(0,0,0,0.4)
    if (S == 1)
      stop("unsupported model with ", sQuote("S=1"))

    r <- polyroot(c(1, rbind(array(0, dim=c(S-1, x$arma[3])), sarcoefs)))

    tmp <- findInterval(abs(Arg(r)), 
      c(rbind(wseas - width[2], wseas + width[2])), rightmost.closed=TRUE)
    tmp <- tmp %% 2
    
    idaux <- which(tmp == 0)
    id.tc <- idaux[which(1/Mod(r[idaux]) > min.modulus)]
    id.trans <- idaux[which(1/Mod(r[idaux]) <= min.modulus)]
    id.seas <- which(tmp != 0)

    rtrend.s <- c(rtrend.s, r[id.tc])
    rtrans <- c(rtrans, r[id.trans])
    rseas.s <- c(rseas.s, r[id.seas])
  }

  # build stationary and total polynomials

  ptrend.s <- roots2poly(rtrend.s)
  ptrend.total <- polyprod(ptrend.ns, ptrend.s)

  #NOTE non-seasonal transitory component is not considered
  ptrans <- roots2poly(rtrans)

  pseas.s <- roots2poly(rseas.s)
  pseas.total <- polyprod(pseas.ns, pseas.s)

  ptotal <- polyprod(c(1, -x$model$Delta), c(1, -x$model$phi))

  structure(list(
    arcoefs = arcoefs, sarcoefs = sarcoefs, d = d, D = D, period = S, width = width,
    roots.nonstationary = list(trend = rtrend.ns, seasonal = rseas.ns),
    roots.stationary = list(trend = 1/rtrend.s, transitory = 1/rtrans, seasonal = 1/rseas.s),
    polys.nonstationary = list(trend = ptrend.ns, seasonal = pseas.ns),
    polys.stationary = list(trend = ptrend.s, transitory = ptrans, seasonal = pseas.s),
    #polys.total = list(trend = ptrend.total, seasonal = pseas.total),
    total = ptotal, 
    trend = ptrend.total, transitory = ptrans, seasonal = pseas.total),
    class="tsdecARroots")
}
