\name{coord.hyp.basis}
\alias{coord.hyp.basis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Internal function to find the basis of a subspace }
\description{
  If a dimension reduction object has \eqn{p} dimension, returns a \eqn{p \times r}
  matrix whose columns space a subspace specified by a formula.
}
\usage{
coord.hyp.basis(object, spec, which = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{dr} object }
  \item{spec}{ A one-sided formula, see below }
  \item{which}{ either +1 or \eqn{-1}{-1}; see below }
}
\details{
The workings of this function is best explained by an example.  Suppose
the \code{dr} object was created with the formula \code{y~x1+x2+x3+x4}, so
we have \eqn{p=4} predictors.  A matrix that spans the subspace of 
\eqn{R^4} specified
by Span(x1,x2,x3,x4) is simply the identity matrix of order 4.

This function will return a subset of the columns of this identity matrix,
as determined by \code{spec}.  For example, if \code{spec = ~.-(x3+x4)}, the 
function returns the columns corresponding to x1 and x2 if \code{which=+1} or the
columns corresponding to x3 and x4 if which=\eqn{-1}{-1}.  Similarly, if 
\code{spec=~x1+x2}, the same matrices will be returned.
}
\value{
  A matrix corresponding to the value of \code{spec} and \code{which} given.
}
\author{Sanford Weisberg, sandy@stat.umn.edu}

\seealso{ \code{\link{dr.coordinate.test}} }
\examples{
data(ais)
s1 <- dr(LBM~log(Ht)+log(Wt)+log(RCC)+log(WCC)+log(Hc)+log(Hg),
         data=ais,method="sir")
coord.hyp.basis(s1,~.-log(Wt)-log(Hg))
 }
            
\keyword{ internal }
