\name{dr.slices}
\alias{dr.slices}
\alias{dr.slice1}
\alias{dr.slice2}
\alias{dr.slice.1d}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Divide a vector into slices of approximately equal size }
\description{
  Divides a vector into slices of approximately equal size.
}
\usage{
dr.slices(y, nslices)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{y}{a vector of length n or an n by p matrix}
  \item{nslices}{the number of slices, no larger than n, or a vector of
  p numbers giving the number of slices in each direction.  If y has p 
  columns and nslices is a number, then the number of slices in each direction
  is the smallest integer greater than the p-th root of nslices.}
}
\details{
  If y is an n-vector, order y.  The goal for the number of observations per slice
  is m, the smallest integer in nslices/n. Allocate the first m observations to
  slice 1.  If there are duplicates in y, keep adding observations to the first
  slice until the next value of y is not equal to the largest value in the 
  first slice.  Allocate the next m values to the next slice, and again check 
  for ties.  Continue until all values are allocated to a slice.  This does not
  guarantee that nslices will be obtained, nor does it guarantee an equal number
  of observations per slice.  This method of choosing slices is invariant under
  rescaling, but not under multiplication by -1, so the slices of y will not
  be the same as the slices of -y.
  
  If y is a matrix of p columns, slice the first column as described above.  Then,
  within each of the slices determined for the first column, slice based on the
  second column, so that each of the "cells" has approximately the same number
  of observations.  Continue through all the columns.  This method is not
  invariant under reordering of the columns, or under multiplication by -1.
}
\value{
  \item{slice.sizes} {gives the number
  of observations per slice.  }
  \item{slice.indicator} {is an n-vector giving the
  slice number of each observation.  }
  \item{nslices} {gives the actual number of slices produced}
}
\references{ This function produces the same slice indicators as does the
equivalent function in Arc, Cook and Weisberg (1999), Applied Regression Including
Computing and Graphics, New York:  Wiley. }
\author{ Sanford Weisberg, <sandy@stat.umn.edu> }

\seealso{  \code{\link{dr}} }

\keyword{ regression }% __ONLY ONE__ keyword per line
