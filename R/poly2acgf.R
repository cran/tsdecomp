
poly2acgf <- function(x, type = c("roots2poly", "acov2ma"), 
  tol = 1e-16, maxiter = 100, init.tol = 1e-05, init.maxiter = 100)
{
  type <- match.arg(type)

  if (type == "roots2poly")
  {
    ##NOTE 
    #based on  algorithm 'pu2ma' in software SSMMATLAB by Gómez, V. 

    r <- polyroot(x)
    tmp <- sqrt(r^2 - 4)
    r <- c((r + tmp)/2, (r - tmp)/2)
    r <- r[order(Mod(r))]
    r <- r[seq_len(length(r)/2)]

    res <- list(coef=roots2poly(1/r), sigma2=x[length(x)]/Re(prod(-r)))

  } else # type == "acov2ma"
  {
    # this method may be less accurate in the presence of non-invertible roots
    # (non-invertible roots are not discarded)

    z <- c(.matrix.poly2acgf[seq_along(x),seq_along(x)] %*% x)
    res <- acov2ma(x=z, tol=tol, maxiter=maxiter, 
      init=acov2ma.init(z, tol=init.tol, maxiter=init.maxiter)$coef)

    r <- 1/polyroot(res$coef)
  }

  structure(
    list(roots=r, coef = res$coef, sigma2 = res$sigma2),
    class="tsdecMAroots")
}

print.tsdecMAroots <- function(x, units = c("radians", "degrees", "pi"), digits = 4, echo = TRUE, ...)
{
  units <- match.arg(units)
  pi2 <- 2*pi

  r <- x$roots

  A <- Arg(r)
  a <- abs(A)
  a[a<1e-08] <- 0
  a[a>1e+06] <- -Inf
  a[a>1e+06] <- Inf	
  id <- which(A < -1e-08)
  if (length(id) > 0)
    a[id] <- pi2 - a[id]
  p <- pi2/a

  if (units == "degrees") {
    a <- a*180/pi
  } else
  if (units == "pi")
  {
##FIXME check
S <- 1
    # rewrite in terms of pi (not the best name, the units are actually radians)
    ref <- round((S/2)*(a/pi), 12)
    id <- which(!(ref %in% seq_len(S)))
    a[id] <- round(a[id], 3)
    id <- which(ref %in% seq_len(S))
    a[id] <- sprintf("%dpi/%d", ref[id], as.integer(S/2))
    id <- which(ref == 1)
    a[id] <- sprintf("pi/%d", as.integer(S/2))
    id <- which(ref == S/2)
    a[id] <- "pi"
  }

  #Cycles.per.Year=S/p
  m <- data.frame(Root=r, Modulus=Mod(r), Argument=a, Period=p)
  #m <- m[order(m[,"Cycles.per.Year"]),]
  m <- m[order(m[,"Argument"]),]
  rownames(m) <- seq_len(nrow(m))

  if (echo)
  {
    cat(paste("Roots of MA polynomial\n----------------------\n"))
    polychar <- polystring(x$coef, varchar="L", brackets=TRUE, ndec=2)
    cat(paste0(polystring(x$coef, varchar="L", brackets=TRUE, ndec=2), 
      "e_t,  e_t ~ IID(0, ", round(x$sigma2, digits=digits), ") \n\n"))

    print(m, digits=digits, ...) #right=FALSE
  }

  invisible(m)
}
