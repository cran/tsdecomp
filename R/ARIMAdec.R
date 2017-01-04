
ARIMAdec <- function(x, mod, width = c(0.035, 0.035), min.modulus = 0.4, 
  extend = 16, drift = FALSE, optim.tol = 1e-04, ...)
{
  #do.call("roots.allocation", )
  p <- roots.allocation(mod, width=width, min.modulus=min.modulus)

  psp <- pseudo.spectrum(mod, p)

  cd <- canonical.decomposition(psp$num$trend, psp$den$trend, 
    psp$num$trans, psp$den$trans, psp$num$seas, psp$den$seas,
    psp$quotient, optim.tol = optim.tol, ...)

  fres <- filtering(x=x, mod=mod, 
    #ma=mod$model$theta[seq_len(mod$arma[2]+mod$arma[4]*mod$arma[5])], 
    trend=list(ar=p$trend, ma=cd$trend$coef, sigma2=cd$trend$sigma2),
    transitory=list(ar=p$trans, ma=cd$trans$coef, sigma2=cd$trans$sigma2),
    seasonal=list(ar=p$seas, ma=cd$seas$coef, sigma2=cd$seas$sigma2),
    irregular.sigma2=cd$irregular.sigma2,
    extend=extend, drift=drift)

  structure(list(ar=p, spectrum=psp, ma=cd, 
    xextended=fres$xext, filters=fres$filters, components=fres$comp), 
    class="ARIMAdec")
}

print.ARIMAdec <- function(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
{
  print(x$ar, units=units)
  cat("\n")
  print(x$ma, units=units)
}

plot.ARIMAdec <- function(x, ...)
{
  #plot(x$components, ...)
  plot(structure(list(components = x$components), class="tsdecFilter"), ...)
}
