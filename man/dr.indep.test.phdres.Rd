\name{dr.indep.test.phdres}
\alias{dr.indep.test.phdres}
\alias{cov.ew.matrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Test of independence for pHd based on residuals }
\description{
  Compute the p-value for the hypothesis that the central
  subspace for OLS residuals has dimension 0 against the alternative that it
  is greater than zero.  This is an
  internal method not usually needed by the user.
}
\usage{
dr.indep.test.phdres(object, stat)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ phdres object }
  \item{stat}{ value of the test statistic }
}

\references{ R. D. Cook (1998) Regression Graphics, New York:  Wiley. }
\author{ Sanford Weisberg <sandy@stat.umn.edu> }

\seealso{ \code{\link{dr}} }
\keyword{regression}
