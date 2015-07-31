\name{dr.permutation.test}
\alias{dr.permutation.test}
\alias{dr.permutation.test.statistic}

\title{Permutation tests of dimension for dr}
\description{
Approximates marginal dimension test significance levels for sir, save,
and phd by sampling from the permutation distribution.}
\usage{
dr.permutation.test(object, npermute=50,numdir=object$numdir)

}

\arguments{
 \item{object}{a dimension reduction regression object created by dr}
 \item{npermute}{number of permutations to compute, default is 50}
 \item{numdir}{maximum permitted value of the dimension, 
with the default from the object}
}

\details{
The method approximates significance levels of the marginal dimension
tests based on a permutation test.  The algorithm:  (1) permutes the
rows of the predictor but not the response; (2) computes marginal
dimension tests for the permuted data; (3) obtains significane levels 
by comparing the observed statsitics to the permutation distribution.

The method is not implemented for ire.
}
\value{
Returns an object of type \sQuote{dr.permutation.test} that can be printed or
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
}

\keyword{regression}%-- one or more ...
