\name{dr.fit}
\alias{dr.fit}
\alias{dr.fit.default}
\title{ Fit dimension reduction regression }
\description{
  Internal generic function that estimates the central subspace.  
}
\usage{
dr.fit(object, numdir=4, ...)

}
\arguments{
  \item{object}{ dimension reduction regression object }
  \item{numdir}{ maximum number of dimensions to consider }
  \item{tol}{ tolerance passed to singular value decomposition }
  \item{\dots}{ other arguments passed to dr.fit.M }
}

\details{
These functions will not
  typically be called directly by the user.  At present, the same dr.fit method works 
  for all dimension reduction methods implemented in this package, but one could
  potentially write a special dr.fit method if needed.
  
  The general outline of this method is as follows.  (1)  A matrix M is computed
  by a call to dr.fit.M(object,...), such that the columns of M are estimated to 
  fall in the subspace of interest (either the central subspace or the central mean
  subspace). (2)  If M is square, its eigenvalues and eigenvectors are computed; if
  M is not square, the eigenvalues of M'M are computed. (3) M was computed with scaled
  and centered predictors.  The eigenvectors are backtransformed to the original
  scale.} 

\value{
  \item{evectors }{ordered eigenvectors that describe the estimates of the 
  dimension reduction subspace}
  \item{evalues }{ordered eigenvalues}
  \item{numdir}{number of eigenvalues}
  \item{raw.evectors}{eigenvectors of the rotated data}
  \item{M}{The kernel matrix.}
}
%\references{ ~put references to the literature/web site here ~ }
\author{Sanford Weisberg <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}} }
\keyword{regression}
\keyword{internal}
