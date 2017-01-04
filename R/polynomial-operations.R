
polystring <- function(x, varchar = "x", brackets = FALSE, ndec = 2, emptychar = "")
{
  if (length(x) > 0)
  {
    x <- round(x, ndec)
    vals <- apply(rbind(
      sapply(sign(x[-1]), FUN=function(x) if (x==1) " + " else " - "), abs(x[-1])), 
      MARGIN=2, paste0, collapse="")
    res <- x[1]
    if (length(vals) > 0) 
      res <- paste0(res, vals[1], varchar, collapse="")
    if (length(vals) > 1) 
      res <- paste0(res, paste0(vals[-1], varchar, "^", seq_along(vals)[-1], collapse=""))
    if (brackets)
      res <- paste0("(", res, ")")
  } else 
  res <- emptychar

  gsub(sprintf("\\s1%s", varchar), sprintf(" %s", varchar), res)
}

polyeval <- function(p, x)
  if (length(p) == 1) p else sum(p * c(1, x^seq_along(p[-1])))

polyprod <- function(x, y, tol = 1.490116e-08) 
{
  #based on stats::convolve(x, rev(y), type="open")

  #tol <- sqrt(.Machine$double.eps)
  n <- length(x)
  x <- c(rep.int(0, length(y) - 1), x)
  n <- length(y <- c(rev(y), rep(0, n - 1)))
  x <- Re(stats::fft(stats::fft(x) * Conj(stats::fft(y)), inverse=TRUE))/n

  if (tol != 0)
    x[abs(x) < tol] <- 0

  x
}

polydiv <- function(x, y)
{
  # based on Pollock (1999) algorithm (4.46)

  nx <- length(x) - 1
  ny <- length(y) - 1

  if (nx < ny)
    return(list(quotient = 0, remainder = x))

  quotient <- c(rep(0, nx-ny), x[nx+1] / y[ny+1])

  for (j in seq.int(nx-1, by=-1, length.out=(nx-ny)))
  {
    s <- min(j, nx-ny) + 1
    ref <- j - ny + 1
    tmp <- sum(quotient[seq.int(ref+1, s)] * y[seq.int(ny, by=-1, length.out=s-ref)])
    quotient[ref] <- (x[j+1] - tmp) / y[ny+1]
  } 

  remainder <- x - polyprod(y, quotient)
  remainder <- remainder[seq_len(nx+1-length(quotient))]
  if (length(remainder) == 0)
    remainder <- 0

  list(quotient = quotient, remainder = remainder)
}

roots2poly <- function(x)
{
  #based on polynom::poly.from.zeros()

  p <- 1
  for (r in x)
    p <- c(p, 0) - c(0, p/r)
  p <- Re(p)

  p
}
