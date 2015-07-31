\name{dr.pvalue}
\alias{dr.pvalue}
\alias{wood.pvalue}
\alias{bentlerxie.pvalue}

\title{ Compute the Chi-square approximations to a weighted
sum of Chi-square(1) random variables. }
\description{
  Returns an approximate quantile for a weighted sum of independent
  \eqn{\chi^2(1)} random variables.
}
\usage{
dr.pvalue(coef,f,chi2approx=c("bx","wood"),...)

bentlerxie.pvalue(coef, f)

wood.pvalue(coef, f, tol=0.0, print=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{ a vector of nonnegative weights }
  \item{f}{ Observed value of the statistic }
  \item{chi2approx}{Which approximation should be used?}
  \item{tol}{tolerance for Wood's method.}
  \item{print}{Printed output for Wood's method}
  \item{...}{Arguments passed from \code{dr.pvalue} to 
wood.pvalue.}
}
\details{
For Bentler-Xie, we approximate \eqn{f} by \eqn{c \chi^2(d)} for values of \eqn{c}
and \eqn{d} computed by the function.  The Wood approximation is more 
complicated.}
\value{
Returns a data.frame with four named components:
\item{test}{The input argument \code{f}.}
\item{test.adj}{For Bentler-Xie, returns \eqn{cf}; for Wood, returns \code{NA}.}
\item{df.adj}{For Bentler-Xie, returns \eqn{d}; for Wood, returns \code{NA}.}
\item{pval.adj}{Approximate p.value.}
}
\references{ Peter M. Bentler and Jun Xie (2000), Corrections to test
statistics in principal Hessian directions.  \emph{Statistics and
Probability Letters}, 47, 381-389. 

Wood, Andrew T. A. (1989)
An \eqn{F} approximation to the distribution of a linear combination of 
chi-squared variables.
\emph{Communications in Statistics: Simulation and Computation}, 18, 1439-1456.}
\author{ Sanford Weisberg <sandy@stat.umn.edu> }


\keyword{ regression }% at least one, from doc/KEYWORDS
