
ARMAacov <- function(ar = numeric(0), ma = numeric(0), lag.max = max(p, q + 1), 
  sigma2 = 1)
{
  ##NOTE based on stats::ARMAacf

  p <- length(ar)
  q <- length(ma)

  if (!p && !q)
    stop("empty model supplied")
  r <- max(p, q + 1)
  #r <- lag.max
  if(p > 0)
  {
    #if(r > 1)
    #{
      if(r > p) { ## pad with zeros so p >= q+1
        ar <- c(ar, rep(0, r - p))
        p <- r
      }
      A <- matrix(c(rep(c(1, -ar, rep(0, p+1)), p), 1, -ar), nrow=p+1, byrow=TRUE)
      A[,seq_len(p)] <- A[,seq_len(p)] + A[,(2L * p + 1L):(p + 2L)]
#debug
#original code in stats::ARMAacf
#a <- matrix(0, p + 1L, 2L * p + 1L)
#ind <- as.matrix(expand.grid(1L:(p + 1), 1L:(p+1)))[, 2L:1L]
#ind[,2] <- ind[,1L] + ind[,2L] - 1L
#a[ind] <- c(1, -ar)
#a[,1L:p] <- a[,1L:p] + a[,(2L * p + 1L):(p + 2L)]
#stopifnot(all.equal(a, A))
      rhs <- c(1, rep(0, p))
      if(q > 0) 
      {
        psi <- c(1, stats::ARMAtoMA(ar, ma, q))
        theta <- c(1, ma, rep(0, q+1L))
        for(k in 1 + 0:q) #seq_len(q)
          rhs[k] <- sum(psi * theta[k + 0:q])
      }
      rhs <- rhs * sigma2
      ind <- (p+1):1
      Acf <- solve(A[ind, ind], rhs)
      #Acf <- Acf[-1L]/Acf[1L]
      #Acf <- Acf[-1L]
    #} else { #Acf <- ar
    #}

    if(lag.max > p)
    {
      Acf <- c(Acf, stats::filter(rep(0, lag.max - p), ar, "recursive", init = rev(Acf[-1L])))
    }
    #Acf <- c(1, Acf[1L:lag.max])
  } else 
  if(q > 0) {
    x <- c(1, ma)
    Acf <- stats::filter(c(x * sigma2, rep(0, q)), rev(x), sides=1)[-seq_len(q)]
    if(lag.max > q)
      Acf <- c(Acf, rep(0, lag.max - q))
    #Acf <- Acf/Acf[1L]
  }

  names(Acf) <- paste0("lag", 0:lag.max)
  Acf
}
