
acov2ma.init <- function(x, tol = 0.00001, maxiter = 100)
{
  # Pollock (1999) procedure (17.35)

  n <- length(x) - 1

##NOTE failure to converge may be due to non-invertible polynomial (change roots)

  sigma2 <- x[1]
  mus <- rep(0, n)

  M <- matrix(seq_len((n-1)*(n-1)), nrow = n-1, ncol = n-1)
  id1 <- M[n-row(M)>=col(M)]
  id2 <- unlist(mapply(seq.int(2, n), FUN=function(i) seq.int(i, n)))
  M[] <- 0
  g2n <- x[seq.int(2,n)]
  gnp1 <-  x[n+1]

  iter <- 0
  conv <- FALSE
  while (!conv && iter < maxiter)
  {
    sigma2.old <- sigma2
    M[id1] <- mus[id2]
    mus <- c(g2n/sigma2 - M %*% mus[-n], gnp1/sigma2)
    sigma2 <- x[1] / (1 + sum(mus^2))

    iter <- iter + 1
    conv <- abs((sigma2.old-sigma2)/sigma2) < tol
  }

  list(convergence = conv, iter=iter, 
    coef = c(1, mus), sigma2 = x[1] / (1 + sum(mus^2)))
}

acov2ma <- function(x, tol = 1e-16, maxiter = 100, init = NULL)
{
  # Pollock (1999) procedure (17.39)

  n <- length(x)

  # initial values of the parameters

  if (is.null(init)) {
    thetas <- c(sqrt(x[1]), rep(0, n-1))
  } else
    thetas <- init

  MH <- Mp <- matrix(0, nrow = n, ncol = n)
  A <- row(MH) + col(MH) - 1
  MH.id1 <- which(A <= n)
  MH.id2 <- A[MH.id1]
  Mp.id1 <- which(upper.tri(Mp, diag = TRUE))
  Mp.id2 <- stats::toeplitz(seq_len(n))[upper.tri(Mp, diag = TRUE)]

  conv <- FALSE
  iter <- 0

  while (!conv && iter < maxiter)
  {
    # build the matrix of derivatives
    # (sigma^2 normalized to 1)

    MH[MH.id1] <- thetas[MH.id2]
    Mp[Mp.id1] <- thetas[Mp.id2]
    MDeriv <- MH + Mp

    # evaluate the function (sigma^2 normalized to 1)

    fvals <- x
    for (i in seq_len(n))
    {
      fvals[i] <- fvals[i] - 
        sum(thetas[seq.int(i, n)] * thetas[seq_len(n-i+1)])
    }

    # find the updating vector

    steps <- solve(MDeriv, fvals)

    # update parameters

    thetas <- thetas + steps

    # check for convergence

    stepsNorm <- crossprod(steps) # sum(steps^2)
    thetasNorm <- crossprod(thetas) #sum(thetas^2)
    #convergence <- if ((stepsNorm/thetasNorm) > tol) FALSE else TRUE
    conv <- (stepsNorm/thetasNorm) < tol
    iter <- iter + 1

  } # end while loop

  # renormalise (theta_0 = 1)

  thetas <- thetas / thetas[1]

  list(convergence = conv, iter = iter, coef = thetas, 
    sigma2 = x[1] / (1 + sum(thetas[-1]^2)))
}
