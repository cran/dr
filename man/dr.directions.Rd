\name{dr.directions}
\alias{dr.direction}
\alias{dr.directions}
\alias{dr.direction.default}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Directions selected by dimension reduction regressiosn }
\description{
  Dimension reduction regression returns a set of p orthogonal direction
  vectors each of length p, the first d of which are estimates a basis of a 
  d dimensional central subspace.  The function returns the estimated directions 
  in the original n dimensional space for plotting.}
}
\usage{
dr.direction(object, which, norm, x)
dr.directions(object, which, norm, x)
dr.direction.default(object, which=1:object$numdir,norm=F,x=dr.x(object))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ a dimension reduction regression object }
  \item{which}{ select the directions wanted, default is all directions }
  \item{norm}{if TRUE, direction vectors are normalized to length 1, otherwise
  their length is arbitrary}
  \item{x}{select the X matrix, the default is dr.x(object)}
  \item{...}{additional arguments are passed to dr.direction.default}
}
\details{
  Dimension reduction regression is used to estimate a basis of the central
  subspace of a regression.  If there are p predictors, the dimension
  reduction regression object includes a p by p matrix of C of eigenvectors.
  This method returns (X-m1')C where m is the vector of column means of X.  If
  X is equal to the original matrix of predictors given by dr.x(object), then
  this gives the directions in the coordinates of the orginal n dimensional
  space.  These directions are used in graphical methods and elsewhere.
}
\value{
  Returns a matrix.  The same function has two names.
}
\references{ See R. D. Cook (1998).  Regression Graphics.  New York:  Wiley. }
\author{Sanford Weisberg <sandy@stat.umn.edu> }


\seealso{  \code{\link{dr}} }

\examples{
library(dr)
data(ais)
attach(ais)  # the Australian athletes data
#fit dimension reduction using sir
m1 <- dr(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8)
summary(m1)
dr.directions(m1)
}
\keyword{ regression }% __ONLY ONE__ keyword per line
