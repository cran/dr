\name{dr.permutation.test}
\alias{dr.permutation.test}
\alias{dr.permutation.test.statistic}
\alias{dr.permutation.test.statistic.default}
\alias{dr.permutation.test.statistic.phd}
\alias{dr.permutation.test.statistic.phdres}
\alias{dr.permutation.test.statistic.phdy}
\alias{print.dr.permutation.test}
\alias{summary.dr.permutation.test}

\title{ Dimension Reduction Regression Functions}
\description{
These functions require a dimension reduction regression object as input to
produce output of various types.
}
\usage{
dr.permutation.test(object, npermute=50,numdir=object$numdir,permute.weights=TRUE)

}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{object}{a dimension reduction regression object created by dr}
 \item{npermute}{number of permutations to compute, default is 50}
 \item{numdir}{maximum permitted value of the dimension, with the default from
the object}
 \item{permute.weights}{If TRUE, permute weights as well as data in the
 permutation test.  If FALSE, do not permute the weights.}
}
\details{
Approximates significance levels for tests of dimension using a permutation test.  
}
\value{
Returns an object of type 'dr.permutation.test' that can be printed or
summarized to give the summary of the test.
}
\references{See www.stat.umn.edu/arc/addons.html, and then select the article
on dimension reduction regression or inverse regression.}
\author{ Sanford Weisberg, sandy@stat.umn.edu}


\seealso{ \code{\link{dr}} }

\examples{
data(ais)
attach(ais)  # the Australian athletes data
#fit dimension reduction regression using sir
m1 <- dr(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8)
summary(m1)
dr.permutation.test(m1,npermute=100)
plot(m1)
dr.coplot(m1)
}

\keyword{regression}%-- one or more ...
