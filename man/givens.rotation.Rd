\name{givens.rotation}
\alias{givens.rotation}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Create givens rotation matrix}
\description{
For a given angle theta, returns a p by p Givens rotation matrix.
}
\usage{
givens.rotation(theta, p=2, which=c(1, 2))
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{theta}{an angle in radians}
 \item{p}{the dimension of the matrix to be produced}
 \item{which}{two numbers between 1 and p giving the columns/rows for the
nonzero elements of the result.}
}
\value{
Returns a p by p matrix z of all zeroes, except z[which,which] has elements
cos(theta), -sin(theta), sin(theta) and cos(theta), in column order.
}
\references{
Gene H. Golub and Charles F. Van Loan (1989).  Matrix Computations, Second
Edition.  Baltimore:  Johns Hopkins Press, p. 202.}
\author{sandy@stat.umn.edu}

\examples{
 givens.rotation(1,4,c(1,3))
}
\keyword{ regression }%-- one or more ...
