\name{dr.y}
\alias{dr.y}
\alias{dr.y.default}
\alias{dr.y.phd}
\alias{dr.y.phdres}
\alias{dr.y.phdy}
\alias{dr.y.phdq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Generic function to return the response for a 'dr' fit.
}
\usage{
dr.y(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A 'dr' object}
}
\details{
 This returns the response variable for a 'dr' fit.  For most, dr methods, this
 is just the response as specified in the call to 'dr.'  For phd and phdres, this
 returns the residuals from the OLS regression of the original response on the
 design matrix.  For 'phdy', it returns the response minus its mean.
}
\value{
  A vector.
}
\author{Sanford Weisberg, sandy@stat.umn.edu }

\seealso{ ~~objects to See Also as \code{\link{dr}}. }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(object) {UseMethod("dr.y")}
}
\keyword{ internal}% at least one, from doc/KEYWORDS
