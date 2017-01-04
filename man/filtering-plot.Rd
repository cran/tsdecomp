\name{plot.tsdecFilter}
\alias{plot.tsdecFilter}

\title{Plot Method for \code{tsdecFilter} Objects}

\description{
Plot the time series containing the components in a \code{tsdecFilter} object.
}

\usage{
\method{plot}{tsdecFilter}(x, select = colnames(X), 
  overlap.trend = FALSE, args.trend = list(col = "black"),
  set.pars = list(mar = c(0, 3, 0, 3), oma = c(4, 0, 2, 0), mfrow = c(nplot, 1)),
  main = NULL, range.bars = TRUE, ..., 
  col.range = "light gray", 
  args.xlab = list(text = "time", side = 1, line = 2), 
  args.ylab = list(side = 3, adj = 0, line = -1),
  xaxis.line = -0.5)
}

\arguments{
\item{x}{an object of class \code{tsdecFilter} returned by \code{\link[tsdecomp]{filtering}}.}
\item{select}{character vector with the labels of the series to be plot. Allowed values are 
\code{c("observed", "trend", "transitory", "seasonal", "sadj", "irregular")}.},
\item{overlap.trend}{logical, if \code{TRUE} the trend component is plot over the observed data; 
otherwise, the trend is plot separately.}
\item{args.trend}{a list containing the arguments to be passed to \code{\link[graphics]{lines}}. 
If \code{overlap.trend=TRUE}, these options are used to plot the trend; otherwise, it is ignored.}
\item{set.pars}{settings for \code{\link{par}(.)} when setting up the plot.}
\item{main}{plot main title.}
\item{range.bars}{logical indicating if each plot should have a bar at
its right side which are of equal heights in user coordinates. 
The same as in \code{\link[stats]{plot.stl}}.}
\item{\dots}{further arguments passed to \code{\link[graphics]{plot}}.}
\item{col.range}{colour to be used for the range bars, if plotted.
Note this appears after \code{\dots} and so cannot be abbreviated.}
\item{args.xlab}{arguments to be passed to \code{\link[graphics]{mtext}} 
when setting the title for the x axis.}
\item{args.ylab}{arguments to be passed to \code{\link[graphics]{mtext}} 
when setting the title for the y axis.}
\item{xaxis.line}{the number of lines into the margin at which 
the x axis line will be drawn.}
}

\details{
This function is based on \code{\link[stats]{plot.stl}}.
}

%\value{}

\seealso{
\code{\link[tsdecomp]{filtering}}.
}

%\value{}

%\examples{}

\keyword{hplot}
\keyword{ts}
