\name{dr.x}
\alias{dr.x}
\alias{dr.y}
\alias{dr.z}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Inverse regression term matrix}
\description{
dr.x returns the matrix of data constructed from the formula for a
dimension reduction regression.  dr.y returns the response.
}
\usage{
dr.x(object)
dr.y(object)
dr.z(x,weights=NULL,center=TRUE,rotate=TRUE,decomp="svd")
}
\arguments{
 \item{object}{An dimension reduction regression object}
 \item{x}{a matrix to be centered and scaled OR a dimension reduction object
 from which the model matrix will be extracted}
 \item{weights}{a vector of weights, or NULL if unweighted}
 \item{center}{if TRUE, center the output}
 \item{rotate}{if TRUE, multiply by the inverse square root of the sample
 covariance matrix.}
 \item{decomp}{Decomposition to be used in computing the rotation; the
 default is "svd".}
}
\value{
dr.x returns an \eqn{n \times p}{n by p} matrix of terms, excluding the 
intercept, constructed from the dimension reduction regression object.  dr.y returns the
response. dr.z returns a possibly centered and scaled version of x.
}
\author{Sanford Weisberg, <sandy@stat.umn.edu>}

\seealso{dr}

\examples{
data(ais)
attach(ais)
m1 <- dr(LBM~Ht+Wt+RCC+WCC)
dr.x(m1)
dr.y(m1)
}
\keyword{regression}%-- one or more ...
