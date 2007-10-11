\name{dr.slices}
\alias{dr.slices}
\alias{dr.slices.arc}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Divide a vector into slices of approximately equal size }
\description{
  Divides a vector into slices of approximately equal size.
}
\usage{
dr.slices(y, nslices)

dr.slices.arc(y, nslices)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{a vector of length \eqn{n} or an \eqn{n \times p} matrix}
  \item{nslices}{the number of slices, no larger than \eqn{n}, or a vector of
  \eqn{p} numbers giving the number of slices in each direction.  If \eqn{y} 
  has \eqn{p} 
  columns and nslices is a number, then the number of slices in each direction
  is the smallest integer greater than the p-th root of nslices.}
}
\details{
  If \eqn{y} is an n-vector, order \eqn{y}.  The goal for the number of observations per slice
  is \eqn{m}, the smallest integer in nslices/n. Allocate the first \eqn{m} observations to
  slice 1.  If there are duplicates in \eqn{y}, keep adding observations to the first
  slice until the next value of \eqn{y} is not equal to the largest value in the 
  first slice.  Allocate the next \eqn{m} values to the next slice, and again check 
  for ties.  Continue until all values are allocated to a slice.  This does not
  guarantee that nslices will be obtained, nor does it guarantee an equal number
  of observations per slice.  This method of choosing slices is invariant under
  rescaling, but not under multiplication by \eqn{-1}, so the slices of \eqn{y} will not
  be the same as the slices of \eqn{-y}.  This function was rewritten for Version 2.0.4 of
  this package, and will no longer give exactly the same results as the program Arc.  If you
  want to duplicate Arc, use the function \code{dr.slice.arc}, as illustrated in the 
  example below.

  If \eqn{y} is a matrix of p columns, slice the first column as described above.  Then,
  within each of the slices determined for the first column, slice based on the
  second column, so that each of the \dQuote{cells} has approximately the same number
  of observations.  Continue through all the columns.  This method is not
  invariant under reordering of the columns, or under multiplication by \eqn{-1}.
}
\value{
Returns a named list with three elements as follows:
  \item{slice.indicator}{ordered eigenvectors that describe the estimates of the 
  dimension reduction subspace}
  \item{nslices}{Gives the actual number of slices produced, which may be smaller 
than the number requested.}
  \item{slice.sizes}{The number of observations in each slice.}
}


\references{R. D. Cook and S. Weisberg (1999), \emph{Applied Regression Including
Computing and Graphics}, New York:  Wiley. }
\author{ Sanford Weisberg, <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}} }

\examples{ 
data(ais)
summary(s1 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
                 log(Hc)+log(Ferr), data=ais,method="sir",nslices=8))
# To make this idential to ARC, need to modify slices to match.
summary(s2 <- update(s1,slice.info=dr.slices.arc(ais$LBM,8)))
}
\keyword{ regression }%
