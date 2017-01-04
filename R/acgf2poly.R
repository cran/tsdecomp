
acgf2poly <- function(x)
{
  n <- length(x)

  if (n > 4)
  {
    res <- x
    n2 <- n - 2
    j <- -1
    tmp <- seq.int(2, n-1)
    
    #with n=50:w Warning: integer overflow in 'cumsum'; use 'cumsum(as.numeric(.))' 
    if (n > 40)
      tmp <- as.numeric(tmp)

    for (i in seq.int(3, n-2, 2))
    {
      res[seq_len(n+1-i)] <- res[seq_len(n+1-i)] + j * tmp * x[seq.int(i, n)]
      n2 <- n2 - 2
      j <- -1 * j      
      tmp <- cumsum(tmp[seq_len(n2)])
    }
    if (n%%2 == 0) {
      res[c(1,2)] <- res[c(1,2)] + j * tmp * x[c(n-1,n)]
    } else
      res[1] <- res[1] + j * tmp * x[n]

  } else 
  if (n == 4)
  {
    res <- c(x[1]-2*x[3], x[2]-3*x[4], x[3], x[4])
  } else 
  if (n == 3)
  {
    res <- c(x[1]-2*x[3], x[2], x[3])
  } else { # n <= 2
    res <- x
  }

  res
}
