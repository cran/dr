\name{dr.estimate.weights}
\alias{dr.estimate.weights}
\alias{robust.center.scale}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute estimated weighting toward normality}
\description{
  This function estimate weights to apply to the rows of a data matrix to
  make the resulting weighted matrix as close to multivariate normality as
  possible.  This method is usually not called directly by the user.
}
\usage{
dr.estimate.weights(object, sigma=1, covmethod="mve", nsamples=NULL, ...)
robust.center.scale(x, method,... )
}
\arguments{
  \item{object}{a dimension reduction regression object name, or an n by p matrix}
  \item{sigma}{A tuning parameter, with default 1, usually in the range .2
  to 1.0}
  \item{covmethod}{covmethod is passed as the argument \code{method} 
  to the function \code{cov.rob} in the required package
  \code{lqs}.  The choices are \code{"classical"},
  \code{"mve"} and \code{"mcd"}.  This probably will not work with Splus.
  If classical is selected, the usual estimate of the covariance matrix is
  used, but the center is the medians, not the means.}
  \item{nsamples}{The weights are determined by random sampling from a
  data-determined normal distribution.  This controls the number of samples}
  \item{x}{An \eqn{n \times p}{n by p} data matrix with no missing values}
  \item{method}{see covmethod above}
  \item{\dots}{Additional args are passed to \code{cov.rob}}
}
\details{
  The basic outline is:  (1) Estimate a mean m and covariance matrix S using a
  possibly robust method; (2) For each iteration, obtain a random vector
  from N(m,sigma*S).  Add 1 to a counter for observation i if the i-th row
  of the data matrix is closest to the random vector; (3) return as weights
  the sample faction allocated to each observation.  If you set the keyword
  \code{weights.only} to \code{T} on the call to \code{dr}, then only the
  list of weights will be returned.}
}
\value{
  Returns a list of n weights, some of which may be zero.  
}
\references{ R. D. Cook and C. Nachtsheim (1994), Reweighting to achieve
elliptically contoured predictors in regression.  Journal of the American
Statistical Association, 89, 592--599.}
\author{Sanford Weisberg, sandy@stat.umn.edu}

\seealso{\code{\link{dr.weights}}, \code{\link{lqs}} }

%\examples{
%}
\keyword{ regression}%-- one or more ...
