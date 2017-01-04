
partial.fraction <- function(numerator, den.trend, den.transitory, den.seasonal)
{
  # NOTE this function assumes that the total denominator in the left-hand-side 
  # of the pseudo-spectrum is the product of the three denominators in the 
  # right-hand-side ("den.trend", "den.transitory", "den.seasonal")

  #NOTE length(den.trend) == 1 && den.trend==1 is slightly faster than "identical"
  #"den.trend==1" may be not needed;
  #alternatively pseudo.spectrum could return NULL, but there it is more convenient 
  #to define the polynomial as 1 when the component is not present,
  #changing there the output would require some extra "if" statements

  ref <- paste0(as.numeric(c(length(den.trend) == 1 && den.trend==1, 
    length(den.transitory) == 1 && den.transitory==1, 
    length(den.seasonal) == 1 && den.seasonal==1)), collapse="")

##FIXME do for other combinations
  switch(ref,
  "010" = 
  {
    p1 <- den.seasonal
    p3 <- den.trend
  },
  "000" = 
  {
    p1 <- polyprod(den.transitory, den.seasonal)
    p2 <- polyprod(den.trend, den.seasonal)
    p3 <- polyprod(den.trend, den.transitory)
  },
  "011" = 
  {
    # p1 <- 1
    # there is only one component, 
    # the LHS numerator will be equal to the seeked RHS numerator
    p1 <- polyprod(den.transitory, den.seasonal)
  },
  "101" =
  {
    #p2 <- 1
    p2 <- polyprod(den.trend, den.seasonal)
  },
  "001" =
  {
    p1 <- den.transitory
    p2 <- den.trend
  })

  nam1 <- length(den.trend) - 2
  nbm1 <- length(den.transitory) - 2
  ncm1 <- length(den.seasonal) - 2

  nr <- max(0, nam1+1) + max(0, nbm1+1) + max(0, ncm1+1)

  if (nam1 < 0) {
    m1 <- matrix(nrow=nr, ncol=0)
  } else
  if (nam1 == 0) {
    m1 <- cbind(p1)
  } else
    m1 <- matrix(c(rep(c(p1, rep(0, nr-length(p1)+1)), nam1), p1), nrow=nr)

  if (nbm1 < 0) {
    m2 <- matrix(nrow=nr, ncol=0)
  } else
  if (nbm1 == 0) {
    m2 <- cbind(p2)
  } else
    m2 <- matrix(c(rep(c(p2, rep(0, nr-length(p2)+1)), nbm1), p2), nrow=nr)

  if (ncm1 < 0) {
    m3 <- matrix(nrow=nr, ncol=0)
  } else
  if (ncm1 == 0) {
    m3 <- cbind(p3)
  } else
    m3 <- matrix(c(rep(c(p3, rep(0, nr-length(p3)+1)), ncm1), p3), nrow=nr)

  lhs <- cbind(m1, m2, m3)
  rhs <- c(numerator, rep(0, nr-length(numerator)))
  res <- as.vector(solve(lhs, rhs))

  list(lhs = lhs, rhs = rhs, 
    num.trend = if (ncol(m1) == 0) NULL else res[seq_len(ncol(m1))], 
    num.transitory = if (ncol(m2) == 0) NULL else 
      res[seq.int(ncol(m1)+1, length.out=ncol(m2))],
    num.seasonal = if (ncol(m3) == 0) NULL else 
      res[seq.int(ncol(m1)+ncol(m2)+1, length.out=ncol(m3))])
}
