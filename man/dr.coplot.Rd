\name{dr.coplot}
\alias{dr.coplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Visualize dimension reduction regression using a coplot }
\description{
  This function produces a coplot of the response versus specified
  directions in dimension reduction regression.
}
\usage{
dr.coplot(object, which=1:object$numdir, mark.by.y=F, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a dr object }
  \item{which}{ specify the indices of the directions to be plotted }
  \item{mark.by.y}{ if FALSE, plot the response versus the first specified
  direction conditioning on the second specified direction.  If TRUE, plot
  first direction versus second direction conditioning on the third
  direction, using the response as a marking variable to color the points.}
  \item{\dots}{ other arguments passed to coplot }
}
\details{
The function dr.permutation.test.statistic is not to be called by the user,
but it computes the test statistic for each different type of dr object.
}
\value{
  A coplot.
}
\author{ Sanford Weisberg <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}} }

\keyword{ regression }% at least one, from doc/KEYWORDS
