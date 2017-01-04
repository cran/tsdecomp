
plot.tsdecFilter <- function(x, #labels = colnames(X),
  select = colnames(X), overlap.trend = FALSE, args.trend = list(col = "black"),
  set.pars = list(mar = c(0, 3, 0, 3), oma = c(4, 0, 2, 0), mfrow = c(nplot, 1)),
  main = NULL, range.bars = TRUE, ..., 
  col.range = "light gray", 
  args.xlab = list(text = "time", side = 1, line = 2), 
  args.ylab = list(side = 3, adj = 0, line = -1),
  xaxis.line = -0.5)
{
  ##NOTE this function is based on plot.stl() in package "stats"

  X <- x$components
  if (any(!(select %in% colnames(X))))
    stop("wrong names in argument \'select\'")
  X <- X[,select]
  nplot <- ncol(X) - as.numeric(overlap.trend)

  # fill in argument 'set.pars' with remaining default arguments 
  # that may have been omitted in the call
  #ldef <- eval(formals(plot.tsdecFilter)$set.pars)
  ldef <- list(mar = c(0,3,0,3), oma = c(4,0,2,0), mfrow = c(nplot, 1))
  id <- which(!(names(ldef) %in% names(set.pars)))
  if (length(id))
    set.pars <- c(set.pars, ldef[id])

  ldef <- list(text = "time", side = 1, line = 2)
  id <- which(!(names(ldef) %in% names(args.xlab)))
  if (length(id))
    args.xlab <- c(args.xlab, ldef[id])
  #args.xlab <- c(args.xlab, ldef)[unique(c(names(args.xlab), names(ldef)))]

  ldef <- list(side = 3, adj = 0, line = -1)
  id <- which(!(names(ldef) %in% names(args.ylab)))
  if (length(id))
    args.ylab <- c(args.ylab, ldef[id])

  if (range.bars)
    mx <- min(apply(rx <- apply(X, 2, range), 2, diff))

  grDevices::dev.hold()
  on.exit(grDevices::dev.flush())

  if (length(set.pars))
  {
    oldpar <- do.call("par", as.list(names(set.pars)))
    on.exit(par(oldpar), add = TRUE)
    do.call("par", set.pars)
  }

  j <- 0
  for (i in seq_len(ncol(X)))
  {
    if (overlap.trend && select[i] == "trend") {
        next
    } else {
      plot(X[,i], #type = if (i < nplot) "l" else "h",
        type = if (select[i] == "irregular") "h" else "l",
        bty = if (j == 0) "o" else "u",
        xlab = "", ylab = "", 
        xaxt="n", yaxt="n", ...)
        #axes=FALSE ignores value of bty and no box is shown

      if (overlap.trend && select[i] == "observed") {
        do.call("lines", c(list(x=X[,"trend"]), args.trend))
        label <- "observed and trend"
      } else {
        label <- select[i]
        if (label == "sadj")
          label <- "seasonally adjusted"
          #label <- "seas. adj."
      }

      j <- j + 1
    }

    if (range.bars)
    {
      dx <- 1/64 * diff(ux <- par("usr")[1L:2])
      y <- mean(rx[,i])
      graphics::rect(ux[2L] - dx, y + mx/2, ux[2L] - 0.4*dx, y - mx/2,
       col = col.range, xpd = TRUE)
    }

    if (i == 1 && !is.null(main))
      graphics::title(main, line = 2, outer = par("oma")[3L] > 0)

    if (select[i] == "irregular")
      abline(h=0)

    right <- j %% 2 == 0
    axis(side=2, labels=FALSE, tcl = 0.25, lwd = 0, lwd.ticks = 1)
    axis(side=2, labels=!right, lwd = 0, lwd.ticks = 0, line = -0.5)
    axis(side=4, labels=FALSE, tcl = 0.25, lwd = 0, lwd.ticks = 1)
    axis(side=4, labels=right, lwd = 0, lwd.ticks = 0, line = -0.5)

    tmp <- c(list(text=paste(" ", label)), args.ylab)
    do.call("mtext", tmp)
  }

  axis(side = 1, labels = FALSE, tcl = 0.25, lwd = 0, lwd.ticks = 1)
  axis(side = 1, lwd = 0, lwd.ticks = 0, line = xaxis.line)

  do.call("mtext", args.xlab)

  invisible()
}
