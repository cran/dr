\name{dr.weights}
\alias{dr.weights}
\title{ Estimate weights for elliptical symmetry}
\description{
  This function estimate weights to apply to the rows of a data matrix to
  make the resulting weighted matrix as close to elliptically symmetric
  as possible.  
}
\usage{
dr.weights(formula, data = list(), subset, na.action = na.fail, 
    sigma=1, nsamples=NULL, ...) 

}
\arguments{
  \item{formula}{A one-sided or two-sided formula.  The right hand side is
used to define the design matrix.} 
  \item{data}{An optional data frame.}
  \item{subset}{A list of cases to be used in computing the weights.}
  \item{na.action}{The default is na.fail, to prohibit computations.  
If set to na.omit, the function will return a list of weights of the 
wrong length for use with dr.}
  \item{nsamples}{The weights are determined by random sampling from a
  data-determined normal distribution.  This controls the number of samples.  The
  default is 10 times the number of cases.}
  \item{sigma}{Scale factor, set to one by default; see the paper by
Cook and Nachtsheim for more information on choosing this parameter.}
  \item{\dots}{Arguments are passed to \code{\link{cov.rob}} to compute a
robust estimate of the covariance matrix.}
}

\details{
  The basic outline is:  (1) Estimate a mean m and covariance matrix S using a
  possibly robust method; (2) For each iteration, obtain a random vector
  from N(m,sigma*S).  Add 1 to a counter for observation i if the i-th row
  of the data matrix is closest to the random vector; (3) return as weights
  the sample faction allocated to each observation.  If you set the keyword
  \code{weights.only} to \code{T} on the call to \code{dr}, then only the
  list of weights will be returned.}

\value{
  Returns a list of \eqn{n} weights, some of which may be zero.  
}
\references{ R. D. Cook and C. Nachtsheim (1994), Reweighting to achieve
elliptically contoured predictors in regression.  Journal of the American
Statistical Association, 89, 592--599.}
\author{Sanford Weisberg, sandy@stat.umn.edu}

\seealso{\code{\link{dr}}, \code{\link{cov.rob}} }

\examples{
data(ais)
w1 <- dr.weights(~ Ht +Wt +RCC, data = ais)
m1 <- dr(LBM~Ht+Wt+RCC,data=ais,weights=w1)
}
\keyword{ regression}%-- one or more ...
