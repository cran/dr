\name{dr.x}
\alias{dr.x}
\alias{dr.y}
\alias{dr.z}

\title{ Accessor functions for data in dr objects }
\description{
  Accessor functions for dr objects.
}
\usage{
dr.x(object)
dr.y(object)
dr.z(object) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object that inherits from \code{dr}. }
}

\value{
Returns a component of a dr object.  \code{dr.x} returns the predictor
matrix reduced to full rank by dropping trailing columns; \code{dr.y}
returns the response vector/matrix, and \code{dr.z} returns the centered
and scaled predictor matrix.
}

\author{ Sanford Weisberg, <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}}. }

\keyword{internal}% at least one, from doc/KEYWORDS
