\name{rotplot}
\alias{rotplot}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{draw many 2D projections of a 3D plot}
\description{
This function draws several 2D views of a 3D plot, sort of like a spinning
plot.
}
\usage{
rotplot(x, y, number=9, theta=seq(0, pi/2, length = 9), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{x}{a matrix with 2 columns giving the horizontal axes of the full 3D
plot.}
 \item{y}{the vertical axis of the 3D plot.}
 \item{number}{Number of 2D views of the 3D plot.}
 \item{theta}{a list of rotation angles}
 \item{\dots}{additional arguments passed to coplot}
}
\details{
For each value of theta, draw the plot of cos(theta)*x[,1]+sin(theta)*x[,2]
versus y.
}
\value{
returns a graph object.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Sanford Weisberg, sandy@stat.umn.edu}
%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~
%
%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }
%
\examples{
 data(ais)
 attach(ais)
 m1 <- dr(LBM ~ Ht + Wt + WCC)  % Estimate central subspace via SIR
 rotplot(dr.direction(m1,which=1:2),dr.y(m1),col=markby(Sex))
  }
\keyword{regression}%-- one or more ...
