\name{dr.M}
\alias{dr.M}
\alias{dr.M.ols}
\alias{dr.M.sir}
\alias{dr.M.save}
\alias{dr.M.phd}
\alias{dr.M.phdy}
\alias{dr.M.phdq}
\alias{dr.M.mphdq}
\alias{dr.M.phdres}
\alias{dr.M.msir}
\alias{dr.M.mphd}
\alias{dr.M.mphdres}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  Compute the kernel matrix for a dr method
}
\usage{
dr.M(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of type 'dr'. }
  \item{\dots}{ Additional method-specific arguments; see description below}
}
\details{
  This is a method-specific function to compute the kernel matrix M of a
  dimension reduction method.  For sir, save, msir and save, the additional 
  arguements are 'nslices', the number of slices to be used, and 'slice.info',
  which if present is the a list produced by \code{\link{dr.slices}}.  For ols and
  the phd methods, no additional arguments are required.
  
  New dimension reduction methods can be added to 'dr' by writing a new 'dr.M' 
  method.
}
\value{
Returns a list with components
  \item{M }{The kernel matrix}
  \item{slice.info }{Information about the slices, if used by the dr method.}
  ...
}
\author{ Sanford Weisberg, sandy@stat.umn.edu }

\seealso{ ~~objects to See Also as \code{\link{dr}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(object, ...){UseMethod("dr.M")}
}
\keyword{ internal }% at least one, from doc/KEYWORDS
