\name{wood.pvalue}
\alias{wood.pvalue}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Chi-square approximation to a sum of Chi-squares}
\description{
Returns an approximation to P(coef'X > f) for X=(X1,...,Xk)', a vector of iid
one df chi-squared rvs.  coef is a list of positive coefficients. tol is used
to check for near-zero conditions.
}
\usage{
wood.pvalue(coef, f, tol=0, print=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{coef}{ Vector of positive numbers, usually from the eigenvalues of a
  random matrix}
  \item{f}{ observed value of statistic}
  \item{tol}{  A small number for tolerance}
  \item{print}{ If TRUE, print the output; if FALSE return the pvalue}
}
\details{
  Following Wood (1989), use either a Saterthwaite-Welsh, inverse gamma or
  three-parameter Pearson Type IV approximation, depending on the values of
  the coefficients.
}
\value{
  Returns the tail area probability.
}
\references{
See Wood (1989), Communications in Statistics, Simulation 1439-1456.
Translated from Arc function wood-pvalue.
}

\author{ Sanford Weisberg <sandy@stat.umn.edu>.  Arc version by Dennis
Cook.}

\seealso{ \code{\link{dr}} }

\keyword{ regression }% __ONLY ONE__ keyword per line
