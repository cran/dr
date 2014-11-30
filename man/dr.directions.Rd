\name{dr.directions}
\alias{dr.direction}
\alias{dr.directions}
\alias{dr.direction.default}
\alias{dr.direction.ire}
\alias{dr.basis}
\alias{dr.basis.default}
\alias{dr.basis.ire}

\title{ Directions selected by dimension reduction regressiosn }
\description{
  Dimension reduction regression returns a set of up to \eqn{p} orthogonal direction
  vectors each of length \eqn{p}, the first \eqn{d} of which are estimates a basis of a 
  \eqn{d} dimensional central subspace.  The function returns the estimated directions 
  in the original \eqn{n} dimensional space for plotting.}

\usage{
dr.direction(object, which, x)
dr.directions(object, which, x)
\S3method{dr.direction}{default}(object, which=NULL,x=dr.x(object))

dr.basis(object,numdir)

\method{dr.basis}{ire}(object,numdir=length(object$result))
}

\arguments{
  \item{object}{ a dimension reduction regression object created by dr. }
  \item{which}{ select the directions wanted, default is all directions.
If method is \code{ire}, then the directions depend on the value of the
dimension you select.  If omitted, select all directions. }
  \item{numdir}{The number of basis vectors to return}
  \item{x}{select the X matrix, the default is \code{dr.x(object)}}
}
\details{
Dimension reduction regression is used to estimate a basis of the central
subspace or mean central subspace of a regression.  If there are \eqn{p} 
predictors, the dimension of the central subspace is less than or equal to
\eqn{p}.  These two functions, \code{dr.basis} and \code{dr.direction}, 
return vectors that describe the central subspace in various ways.

Consder \code{dr.basis} first.  If you set \code{numdir=3}, for example, this
method will return a \eqn{p} by 3 matrix whose columns span the estimated
three dimensional central subspace.  For all methods except for \code{ire},
this simply returns the first three columns of \code{object$evectors}.  For
the \code{ire} method, this returns the three vectors determined by a 
three-dimensional solution. Call this matrix \eqn{C}.  The basis is 
determined by back-transforming from centered and scaled predictors to
the scale of the original predictors, and then renormalizing the vectors
to have length one.  These vectors are orthogonal in the inner
product determined by Var(X).

The \code{dr.direction} method return \eqn{XC}, the same space but now a
subspace of the original \eqn{n}-dimensional space.  These vectors are
appropriate for plotting.   
}
\value{
  Both functions return a matrix:  for \code{dr.direction}, the matrix has n rows and
numdir columns, and for \code{dr.basis} it has p rows and numdir columns.
}
\references{ See R. D. Cook (1998).  Regression Graphics.  New York:  Wiley. }
\author{Sanford Weisberg <sandy@stat.umn.edu> }


\seealso{  \code{\link{dr}} }

\examples{
data(ais)
#fit dimension reduction using sir
m1 <- dr(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8, data=ais)
summary(m1)
dr.basis(m1)
dr.directions(m1)
}
\keyword{ regression }% __ONLY ONE__ keyword per line
