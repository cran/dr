\name{drop1.dr}
\alias{drop1.dr}
\alias{dr.step}
\title{ Sequential fitting of coordinate tests using a dr object }
\description{
This function implements backward elimination using a \code{dr} object for which
a \code{dr.coordinate.test} is defined, currently for SIR SAVE, IRE and PIRE.
}
\usage{
dr.step(object,scope=NULL,d=NULL,minsize=2,stop=0,trace=1,...)

\method{drop1}{dr}(object, scope = NULL,  update=TRUE,
test="general",trace=1,\dots)


}
\arguments{
\item{object}{ A \code{dr} object for which \code{dr.coordinate.test} 
is defined, for \code{method} equal to one of \code{sir}, \code{save}
or \code{ire}.}

\item{scope}{A one sided formula specifying predictors that will never
be removed.}

\item{d}{To use \emph{conditional coordinate tests}, specify the dimension
of the central (mean) subspace.  The default is \code{NULL}, meaning no
conditioning. This is currently available only for methods \code{sir}, 
\code{save} without categorical predictors, or for 
\code{ire} with or without categorical predictors.}

\item{minsize}{Minimum subset size, must be greater than or equal to 2.}

\item{stop}{Set stopping criterion:  continue removing variables until
the p-value for the next variable to be removed is less than stop.  The
default is stop = 0.}

\item{update}{If true, the \code{update} method is used to return a
\code{dr} object obtained from \code{object} by updating the formula to
drop the variable with the largest p.value.  This can significantly slow
the computations for IRE but has little effect on SAVE and SIR.}

\item{test}{ Type of test to be used for selecting the next predictor
to remove for \code{method="save"} only. \code{"normal"} assumes normal 
predictors, \code{"general"} assumes elliptically contoured predictors. 
For other methods, this argument is ignored.}

\item{trace}{If positive, print informative output at each step, the default.
If trace is 0 or false, suppress all printing.}

\item{\dots}{Additional arguments passed to \code{dr.coordinate.test}.}
}
\details{
Suppose a \code{dr} object has \eqn{p=a+b} predictors, with \eqn{a} predictors specified in the \code{scope} statement.  
\code{drop1} will compute either marginal coordinate tests (if \code{d=NULL})
or conditional marginal coordinate tests (if \code{d} is positive) for dropping each of the \code{b} predictors not in the scope, and return p.values.
The result is an object created from the original object with the predictor
with the largest p.value removed.

\code{dr.step} will call \code{drop1.dr} repeatedly until 
\eqn{\max(a,d+1)} predictors remain.  }
\value{
As a side effect, 
a data frame of labels, tests, df, and p.values is printed.  If
\code{update=TRUE}, a  \code{dr}
object is returned with the predictor with the largest p.value removed.
}
\references{Cook, R. D. (2004).  Testing predictor contributions in 
sufficient dimension reduction.  \emph{Annals of Statistics}, 32, 1062-1092.

Shao, Y., Cook, R. D. and Weisberg (2007).  Marginal tests with 
sliced average variance estimation.  \emph{Biometrika}.}

\author{Sanford Weisberg, <sandy@stat.umn.edu>, based on the 
\code{drop1} generic function in the
base R. The \code{dr.step} function is also similar to \code{step} in 
base R.}

\seealso{ \code{\link{dr.coordinate.test}}}
\examples{
data(ais)
# To make this idential to ARC, need to modify slices to match by
# using slice.info=dr.slices.arc() rather than nslices=8
summary(s1 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
                 log(Hc)+log(Ferr), data=ais,method="sir",
                 slice.method=dr.slices.arc,nslices=8)) 
# The following will almost duplicate information in Table 5 of Cook (2004).
# Slight differences occur because a different approximation for the
# sum of independent chi-square(1) random variables is used:
ans1 <- drop1(s1)
ans2 <- drop1(s1,d=2)
ans3 <- drop1(s1,d=3)
# remove predictors stepwise until we run out of variables to drop.
dr.step(s1,scope=~log(Wt)+log(Ht))
}
\keyword{regression }% at least one, from doc/KEYWORDS
