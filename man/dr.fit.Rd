\name{dr.fit}
\alias{dr.fit}
\alias{dr.fit.y}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fit dimension reduction regression }
\description{
  Internal function that estimates the central subspace.  This will not
  typically be called directly by the user.  dr.fit.y is the generic name
  of a function that returns the response for the particular dimension
  reduction regression object.  Typically, this just returns a centered
  version of the response, but for phdres it returns OLS residuals.
}
\usage{
dr.fit(object, numdir=4, tol=1e-07, ...)
dr.fit.y(object)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ dimension reduction regression object }
  \item{numdir}{ maximum number of dimensions to consider }
  \item{tol}{ tolerance passed to singular value decomposition }
  \item{\dots}{ other arguments passed to dr.fit.M }
}

\value{
  \item{evectors }{ordered eigenvectors that describe the estimates of the 
  dimension reduction subspace}
  \item{evalues }{ordered eigenvalues}
  \item{numdir}{number of eigenvalues}
  \item{raw.evectors}{eigenvectors of the rotated data}
  \item{decomp}{which decomposition was used?}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Sanford Weisberg <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}} }
\keyword{regression}
