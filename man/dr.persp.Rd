\name{dr.persp}
\alias{dr.persp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Produce a perspective plot of dimension reduction regression objects. }
\description{
  This function uses the sm library to draw a perspective 3D plot of two of
  the directions obtained by dimension reduction regression versus the response.
}
\usage{
dr.persp(object, which=1:2, h=c(0.1, 0.1), ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a dimension reduction regression object }
  \item{which}{ a vector of two directions to use }
  \item{h}{ bandwidths for smoothing passed to sm }
  \item{\dots}{ graphical parameters passed to sm }
}

\value{
  Returns a perspective plot
}
\author{ Sanford Weisberg <sandy@stat.umn.edu>}


\seealso{ \code{\link{dr}}, \code{\link[sm]{sm}} }
\keyword{regression}
