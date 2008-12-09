\name{dr.coordinate.test}
\alias{dr.test}
\alias{dr.coordinate.test}
\alias{dr.joint.test}
\alias{dr.joint.test.default}
\alias{dr.joint.test.ire}
\alias{dr.joint.test.pire}
\alias{dr.coordinate.test.default}
\alias{dr.coordinate.test.sir}
\alias{dr.coordinate.test.save}
\alias{dr.coordinate.test.ire}
\alias{dr.coordinate.test.pire}
\title{Dimension reduction tests}
\description{
Functions to compute various tests concerning the dimension of a
central subspace.
}
\usage{
dr.test(object, numdir, ...)

dr.coordinate.test(object, hypothesis,d,chi2approx,\dots)

\S3method{dr.joint.test}{ire}(object, hypothesis, d = NULL,\dots)
}
\arguments{
\item{object}{ The name of an object returned by a call to \code{dr}.}
\item{hypothesis}{ A specification of the null hypothesis to be tested by
  the coordinate hypothesis.  See details below for options.}
\item{d}{For conditional coordinate hypotheses, specify the dimension of
  the central mean subspace, typically 1, 2 or possibly 3.  If left at the
  default, tests are unconditional.}
\item{numdir}{The maximum dimension to consider.  If not set defaults to 4.}
\item{chi2approx}{Approximation method for p.values of linear combination
  of \eqn{\chi^2(1)}{Chisq(1)} random variables.  Choices are from \code{c("bx","wood")},
  for the Bentler-Xie and Wood approximatations, respectively. The default is either "bx"
  or the value set in the call that created the dr object. }
  \item{...}{Additional arguments.  None are currently available.}
}

\details{
\code{dr.test} returns marginal dimension tests.  
\code{dr.coordinate.test} returns marginal dimension tests (Cook, 2004)
if \code{d=NULL} or conditional dimension tests if \code{d} is a
positive integer giving the assumed dimension of the central
subspace.  The function \code{dr.joint.test} tests the coordinate
hypothesis and dimension simultaneously.  It is defined only for 
ire, and is used to compute the conditional coordinate test.

As an example, suppose we have created a \code{dr} object 
using the formula
\code{y ~ x1 + x2 + x3 + x4}.    
The marginal coordinate hypothesis defined by Cook (2004) tests 
the hypothesis that \code{y} is independent of some of the 
predictors given the other predictors.  For example, one could test 
whether \code{x4} could be dropped from the problem by testing \code{y} 
independent of \code{x4} given \code{x1,x2,x3}.

The hypothesis to be tested is determined by the argument \code{hypothesis}.
The argument \code{hypothesis = ~.-x4} would test the hypothesis of the last
paragraph.  Alternatively, \code{hypothesis = ~x1+x2+x3} would
fit the same hypothesis.

More generally, if \code{H} is a \eqn{p \times q}{p times q} 
rank \eqn{q} matrix, and 
\eqn{P(H)} is the projection
on the column space of \code{H}, then specifying \code{hypothesis = H} will test the 
hypothesis that \eqn{Y} is independent of \eqn{(I-P(H))X | P(H)X}.
}
\value{
Returns a list giving the value of the test statistic and an asymptotic
p.value computed from
the test statistic.  For SIR objects, the p.value is computed in two ways.  The
\emph{general test}, indicated by \code{p.val(Gen)} in the output, assumes only 
that the predictors are linearly related.  The \emph{restricted test}, indicated 
by \code{p.val(Res)} in the output, assumes in addition to the linearity condition
that a constant covariance condition holds; see Cook (2004) for more information
on these assumptions.  In either case, the asymptotic distribution is a linear
combination of Chi-squared random variables.  The function specified by the
\code{chi2approx} approximates this linear combination by a single Chi-squared 
variable.

For SAVE objects, two p.values are also returned.  \code{p.val(Nor)} assumes
predictors are normally distributed, in which case the test statistic is asympotically
Chi-sqaured with the number of df shown.  Assuming general linearly related
predictors we again get an asymptotic linear combination of Chi-squares that
leads to \code{p.val(Gen)}.

For IRE and PIRE, the tests
statistics have an asymptotic \eqn{\chi^2}{Chisq} distribution, so the
value of \code{chi2approx} is not relevant.
}
\references{ Cook, R. D. (2004).  Testing predictor contributions in 
sufficient dimension reduction.  \emph{Annals of Statistics}, 32, 1062-1092.

Cook, R. D. and Ni, L.
(2004).  Sufficient dimension reduction via inverse regression:  A minimum
discrrepancy approach, \emph{Journal of the American Statistical Association},
100, 410-428.

Cook, R. D. and Weisberg, S. (1999). \emph{Applied Regression Including
Computing and Graphics}.  Hoboken NJ:  Wiley.

Shao, Y., Cook, R. D. and Weisberg, S. (2007, in press).  Marginal tests with 
sliced average variance estimation.  \emph{Biometrika}.}

\author{Yongwu Shao for SIR and SAVE and Sanford Weisberg
for all methods, <sandy@stat.umn.edu>}

\seealso{\code{\link{drop1.dr}}, \code{\link{coord.hyp.basis}},
         \code{\link{dr.step}},
         \code{\link{dr.pvalue}}}
\examples{
#  This will match Table 5 in Cook (2004).  
data(ais)
# To make this idential to Arc (Cook and Weisberg, 1999), need to modify slices to match.
summary(s1 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+log(Hc)+log(Ferr),
  data=ais,method="sir",slice.function=dr.slices.arc,nslices=8))
dr.coordinate.test(s1,~.-log(Hg))
#The following nearly reproduces Table 5 in Cook (2004)
drop1(s1,chi2approx="wood",update=FALSE)
drop1(s1,d=2,chi2approx="wood",update=FALSE)
drop1(s1,d=3,chi2approx="wood",update=FALSE)
}
\keyword{regression }% at least one, from doc/KEYWORDS
