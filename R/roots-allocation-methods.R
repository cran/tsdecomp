
print.tsdecARroots <- function(x, units = c("radians", "degrees", "pi"), digits = 4, ...)
{
  units <- match.arg(units)

  pi2 <- 2*pi

  S <- x$period
  d <- x$d
  D <- x$D
  arcoefs <- x$arcoefs
  if (length(arcoefs) > 0) 
    arcoefs <- c(1, arcoefs)
  sarcoefs <- x$sarcoefs
  polychar <- ""

  # differencing filter
  polycharD <- paste0(
    if (d==1) "(1 - L)" else if (d>1) sprintf("(1-L)^%d", d) else "",
    if (D==1) sprintf("(1 - L^%d)", S, D) else if (D>1) sprintf("(1-L^%d)^%d", S, D) else "")

  # AR polynomial character string
  polychar <- polystring(arcoefs, varchar="L", brackets=TRUE, ndec=3)

  # seasonal AR polynomial character string
  if (length(sarcoefs) > 0)
  {
    vals <- apply(rbind(
      sapply(sign(sarcoefs), FUN=function(x) if (x==1) " + " else " - "), round(abs(sarcoefs),2)), 
      MARGIN=2, paste0, collapse="")
    tmp <- paste0("1", vals[1], sprintf("L^%d", S))
    if (length(vals) > 1) 
      tmp <- paste0(tmp, paste0(vals[-1], 
        sprintf("L^%d", seq.int(2*S, S*length(vals), by=S)), collapse=""))
    polychar <- paste0(polychar, paste0("(", tmp, ")"))
  }

  pdecchar <- paste0(
    if (identical(x$trend, 1)) "" else 
      polystring(x$trend, varchar="L", brackets=TRUE, ndec=3),
    if (identical(x$transitory, 1)) "" else  
      polystring(x$transitory, varchar="L", brackets=TRUE, ndec=3),
    if (identical(x$seasonal, 1)) "" else  
      polystring(x$seasonal, varchar="L", brackets=TRUE, ndec=3))
  pdecchar <- gsub(" 1L", " L", pdecchar)

  cat(paste("Roots of AR polynomial\n----------------------\n"))
  cat(paste0(polychar, polycharD, " = ", pdecchar), "\n\n")

  rtrend <- c(x$roots.nonstationary$trend, x$roots.stationary$trend)
  rtrans <- x$roots.stationary$transitory
  rseas <- c(x$roots.nonstationary$seasonal, x$roots.stationary$seasonal)
  labels <- c(rep("trend", length(rtrend)), rep("transitory", length(rtrans)),
    rep("seasonal", length(rseas)))
  if (length(labels) > 0)
  {
    r <- c(rtrend, rtrans, rseas)

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
      # rewrite in terms of pi (not the best name, the units are actually radians)
      ref <- round((S/2)*(a/pi), 12)
      id <- which(!(ref %in% seq_len(S)))
      a[id] <- round(a[id], 3)
      id <- which(ref %in% seq_len(S))
      #if (length(id) > 0)
      a[id] <- sprintf("%dpi/%d", ref[id], as.integer(S/2))
      id <- which(ref == 1)
      a[id] <- sprintf("pi/%d", as.integer(S/2))
      id <- which(ref == S/2)
      a[id] <- "pi"
    }

    m <- data.frame(Component=labels, Root=r, Modulus=Mod(r), Argument=a, 
      Period=p, Cycles.per.Year=S/p)
    m <- m[order(m[,"Cycles.per.Year"]),]
    rownames(m) <- seq_len(nrow(m))

    print(m, digits=digits, ...) #right=FALSE

  } else {
    m <- NULL
    cat(paste("None\n"))
  }

  invisible(m)
}

plot.tsdecARroots <- function(x, xlim, ylim, ...)
{
  ##NOTE this function is based on graphics::pie()

  n1 <- 0 # do not draw a contour line for the unit circle
  n2 <- 100
  eps <- 0

  # complex plane
  # seasonal frequencies
  # based on graphics::pie()

  S <- x$period
  pi2 <- 2*pi

  #S-1 seasonal roots
  if (S == 4) {
    #labels <- c(expression(frac(pi, 2)), expression(pi), expression(frac(3*pi, 2)))
    labels <- c(expression("pi/2"), expression(pi), expression("3pi/2"))
  } else
  if (S == 12)
  {
    labels <- c(expression(frac(pi, 6)), expression(frac(2*pi, 6)),
      expression(frac(3*pi, 6)), expression(frac(4*pi, 6)), expression(frac(5*pi, 6)), 
      expression(pi),
      expression(frac(7*pi, 6)), expression(frac(8*pi, 6)), expression(frac(9*pi, 6)), 
      expression(frac(10*pi, 6)), expression(frac(11*pi, 6)))
  } else
    stop(paste0("unsupported periodicity \"S=", S, "\""))

  rtrend <- c(x$roots.nonstationary$trend, x$roots.stationary$trend)
  rtrans <- x$roots.stationary$transitory
  rseas <- c(x$roots.nonstationary$seasonal, x$roots.stationary$seasonal)

  if (missing(xlim) || missing(ylim))
  {
    r <- c(rtrend, rtrans, rseas)
    xlim <- c(min(-1.15, Re(r)), max(1.15, Re(r)))
    ylim <- c(min(-1.15, Im(r)), max(1.15, Im(r)))
  }

  plot(x=0, xlim=xlim, ylim=ylim, type="n", xlab="", ylab="", ...)

  #based on graphics::pie()
  for (i in seq_len(S-1))
  {
    tt <- seq(pi2*(i/S-x$width[2]), pi2*(i/S+x$width[2]), len=n2)
    P <- list(x=cos(tt), y=sin(tt))
    graphics::polygon(c(P$x, 0), c(P$y, 0), col="gray90", border=NA)

    P <- list(x=cos(pi2*i/S), y=sin(pi2*i/S))
    lines(c(1, 1.05)*P$x, c(1, 1.05)*P$y)
    text(1.1*P$x, 1.1*P$y, labels[i], xpd=TRUE, adj=ifelse(P$x < 0, 1, 0)) #cex=cex.text
  }

  lines(c(1, 1.05), c(0,0))
  text(1.1, 0, "0", xpd=TRUE, adj=0) #cex=cex.text

  # trend
  tt <- seq(0, pi2*x$width[1], len=n2)
  P <- list(x=cos(tt), y=sin(tt))
  graphics::polygon(c(P$x, 0), c(P$y, 0), col="mistyrose", border=NA)
  graphics::segments(0, 0, P$x[length(P$x)], P$y[length(P$y)], col="mistyrose")
  graphics::polygon(c(P$x, 0), c(-P$y, 0), col="mistyrose", border=NA)
  graphics::segments(0, 0, P$x[length(P$x)], -P$y[length(P$y)], col="mistyrose")

  # transitory
  tt <- seq(pi2*(x$width[1]+eps), pi2*(1/S-x$width[2]-eps), len=n2)
  P <- list(x=cos(tt), y=sin(tt))
  graphics::polygon(c(P$x, 0), c(P$y, 0), col="cornsilk", border=NA)

  for (i in seq.int(1, S-1))
  {
    tt <- seq(pi2*(i/S+x$width[2]+eps), pi2*((i+1)/S-x$width[2]-eps), len=n2)    
    P <- list(x=cos(tt), y=sin(tt))
    graphics::polygon(c(P$x, 0), c(P$y, 0), col="cornsilk", border=NA)
  }

  tt <- seq(0, pi2*pi2, len=n1)
  points(cos(tt), sin(tt), pch=".") #col="gray80"

  graphics::segments(-1, 0, 1, 0)
  graphics::segments(0, -1, 0, 1)

  if (length(rtrend) > 0)
  {
    points(Re(rtrend), Im(rtrend), pch=21, col="blue", bg="red", lwd=1, cex=1.2)
    graphics::segments(0, 0, Re(rtrend), Im(rtrend), col="red", lwd=1)
  }

  if (length(rtrans) > 0)
  {  
    points(Re(rtrans), Im(rtrans), pch=21, col="blue", bg="green", lwd=1, cex=1.2)
    graphics::segments(0, 0, Re(rtrans), Im(rtrans), col="green", lwd=1)
  }

  if (length(rseas) > 0)
  {
    points(Re(rseas), Im(rseas), pch=21, col="blue", bg="blue", lwd=1, cex=1.2)
    graphics::segments(0, 0, Re(rseas), Im(rseas), col="blue", lwd=1)
  }
}
