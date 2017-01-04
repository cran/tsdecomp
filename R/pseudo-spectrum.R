
pseudo.spectrum <- function(mod, ar)
{
  ma.total <- c(1, mod$model$theta[seq_len(mod$arma[2]+mod$arma[4]*mod$arma[5])])
#debug
ma <- coef(mod)[seq.int(mod$arma[1]+1, length.out=mod$arma[2])]
sma <- coef(mod)[seq.int(mod$arma[1]+mod$arma[2]+mod$arma[3]+1, length.out=mod$arma[4])]
tmp <- polyprod(c(1, ma), c(1, rbind(array(0, dim=c(mod$arma[5]-1, length(sma))), sma)))
tmp <- c(1, mod$model$theta[seq_len(mod$arma[2]+mod$arma[4]*mod$arma[5])])
stopifnot(all.equal(as.vector(tmp), ma.total))

  if (length(ma.total) > 1) {
    ma.total.bf <- stats::convolve(ma.total, ma.total, type="open")[-seq_along(ma.total[-1])]
  } else
    ma.total <- ma.total.bf <- 1

  den.trend.bf <- stats::convolve(ar$trend, ar$trend, type="open")[-seq_along(ar$trend[-1])]
  den.trans.bf <- stats::convolve(ar$transitory, ar$transitory, type="open")[-seq_along(ar$transitory[-1])]
  den.seas.bf <- stats::convolve(ar$seasonal, ar$seasonal, type="open")[-seq_along(ar$seasonal[-1])]

  # the above coefficients are in terms of 2*cos(wj)
  # transformation into a polynomial in the variable (2*cos(w))

  num.psp.total <- acgf2poly(ma.total.bf)

  #den.psp.trend <- den.trans.trend <- den.psp.seas <- NULL

  if ((ntrend <- max(0, length(den.trend.bf)-1)) > 0) {
    den.psp.trend <- acgf2poly(den.trend.bf)
  } else den.psp.trend <- 1
  
  if ((ntrans <- max(0, length(den.trans.bf)-1)) > 0 ) {
    den.psp.trans <- acgf2poly(den.trans.bf)
  } else den.psp.trans <- 1

  if ((nseas <- max(0, length(den.seas.bf)-1)) > 0 ) {
    den.psp.seas <- acgf2poly(den.seas.bf)
  } else den.psp.seas <- 1

  #total denominator of the pseudo-spectrum
  #this is used only if polynomial division is required
  #"den.psp.total" is not used outside this function, it is obtained and 
  #returned for completeness

  den.psp.total <- polyprod(polyprod(den.psp.trend, den.psp.trans), den.psp.seas)

  if (length(num.psp.total)-1 >= ntrend + ntrans + nseas)
  {
    # division is necessary, otherwhise there would be more elements in 
    # the RHS than equations in the LHS of the system of equations to be solved

    tmp <- polydiv(num.psp.total, den.psp.total)
    quotient <- unname(tmp$quotient)
    num.psp.total <- tmp$remainder
  } else
    quotient <- 0

  # partial fraction decomposition

  pfd <- partial.fraction(num.psp.total, 
    den.psp.trend, den.psp.trans, den.psp.seas)

  structure(list(quotient = quotient, 
    total.numerator = num.psp.total, total.denominator = den.psp.total,
    numerators = list(trend = pfd$num.trend, transitory = pfd$num.transitory, 
    seasonal = pfd$num.seasonal),
    denominators = list(trend = den.psp.trend, transitory = den.psp.trans,
      seasonal = den.psp.seas)), class="tsdecPSP")
}

print.tsdecPSP <- function(x, ...)
{
  nums <- lapply(x$numerators, polystring, ...)
  dens <- lapply(x$denominators, polystring, ...)

  cat("Num/Den = Anum/Aden + Bnum/Bden + Cnum/Cden\n----\n")
  cat("Total numerator (Num):\n")
  cat(polystring(x$total.numerator, ...), "\n")
  cat("Total denominator (Den):\n")
  cat(polystring(x$total.denominator, ...), "\n")

  cat("Trend numerator (Anum):\n")
  cat(nums$trend, "\n")
  cat("Trend denominator (Aden):\n")
  cat(dens$trend, "\n")
  cat("Transitory numerator (Bnum):\n")
  cat(nums$transitory, "\n")  
  cat("Transitory denominator (Bden):\n")
  cat(dens$transitory, "\n")    
  cat("Seasonal numerator (Snum):\n")
  cat(nums$seasonal, "\n")  
  cat("Seasonal denominator (Sden):\n")
  cat(dens$seasonal, "\n")
}
