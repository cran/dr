\name{dr.x}
\alias{dr.x}
\alias{dr.wts}
\alias{dr.qr}
\alias{dr.Q}
\alias{dr.R}
\alias{dr.z}
\alias{dr.y.name}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Accessor functions for dr. objects }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dr.x(object)
dr.wts(object) 
dr.qr(object)
dr.Q(object)
dr.R(object)
dr.z(object) 
dr.y.name(object) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object that inherits from 'dr'. }
}

\value{
 Returns a component of a dr object.  For example, 'dr.qr' returns object\$qr,
 'dr.wts' returns object\$weights, 'dr.Q' returns the Q component of object\$qr,
 'dr.z' returns the scaled and centered design matrix.
}

\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{ ~~objects to See Also as \code{\link{dr}}. }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(object) {object$x}
}
\keyword{internal}% at least one, from doc/KEYWORDS
