\name{dr}
\alias{dr}
\alias{plot.dr}
\alias{print.dr}
\alias{summary.dr}
\alias{print.summary.dr}
\alias{model.matrix.dr}
\alias{dr.weights}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Dimension reduction regression}
\description{
 The function dr implements dimension reduction methods, including SIR, SAVE and pHd.
}
\usage{
dr(formula, data=list(), subset, weights, na.action=na.omit, method="sir", 
     contrasts=NULL,offset=NULL, ...)
dr.weights(formula, data=list(), subset, weights, na.action=na.omit, method="sir", 
     contrasts=NULL,offset=NULL, ...)
 }

\arguments{
 \item{formula}{a symbolic description of the model to be fit. The details of
          the model are the same as for lm. }
 \item{data}{ an optional data frame containing the variables in the model.
           By default the variables are taken from the environment from
           which `dr' is called.}
 \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
 \item{weights}{an optional vector of weights to be used where appropriate.}
 \item{na.action}{a function which indicates what should happen when the data
          contain `NA's.  The default is `na.omit,' which will force
          calculations on a complete subset of cases.}
 \item{method}{This character string specifies the method of fitting.  ``sir"
          specifies sliced inverse regression and ``save" specifies sliced
          average variance estimation.  ``phdy" uses principal hessian
          directions using the response as suggested by Li, and ``phdres" 
          uses the LS residuals as suggested by Cook. Other methods may be
          added}
 \item{contrasts}{an optional list. See the `contrasts.arg' of
          `model.matrix.default'.}
 \item{offset}{Set an offset or NULL}
 \item{\dots}{additional items that may be required or permitted by some 
          methods.
          \code{nslices} is the number of slices used by sir and save.
	  \code{numdir} is the maximum number of directions to compute, with
	  default equal to 4.}
}

\details{
The general regression problem studies \eqn{F(y|x)}, the conditional
distribution of a response \eqn{y} given a set of predictors \eqn{x}.  
This function provides methods for estimating the dimension and central
subspace of a general regression problem.  That is, we want to find a 
\eqn{p \times d}{p by d} matrix \eqn{B} such that 
\deqn{F(y|x)=F(y|B'x)}  
Both the dimension \eqn{d} and the subspace
\eqn{R(B)} are unknown.  These methods make few assumptions.  All the methods
available in this function estimate the unknowns by study of the inverse
problem, \eqn{F(x|y)}.  In each, a kernel matrix \eqn{M} is estimated such
that the column space of \eqn{M} should be close to the central subspace.
Eigenanalysis of \eqn{M} is then used to estimate the central subspace.
Objects created using this function have appropriate print, summary and plot
methods.

Weights can be used, essentially to specify the relative 
frequency of each case in the data.  Empirical weights that make 
the contours of the weighted sample closer to elliptical can be 
computed using \code{\link{dr.weights}}.  
This will usually result in zero weight for some 
cases.  The function will set zero estimated weights to missing.

Several functions are provided that require a dr object as input.  
\code{dr.permutation.tests} uses a permutation test to obtain significance levels
for tests of dimension.  \code{dr.coplot} allows visualizing the results using a
coplot of either two selected directions conditioning on a third and using
color to mark the response, or the resonse versus one direction,
conditioning on a second direction.  \code{plot.dr} provides the default plot
method for dr objects, based on a scatterplot matrix.
}

\value{
dr returns an object that inherits from dr (the name of the type is the
value of the \code{method} argument), with attributes:
  \item{M}{A matrix that depends on the method of computing.  The column space
of M should be close to the central subspace.}
  \item{evalues}{The eigenvalues of M (or squared singular values if M is not
symmetric).}
  \item{evectors}{The eigenvectors of M (or of M'M if M is not square and
symmetric) ordered according to the eigenvalues.}
  \item{numdir}{The maximum number of directions to be found.  The output
value of numdir may be smaller than the input value.}
   \item{ols.coef}{Estimated regression coefficients, excluding the intercept,
for the (weighted) LS fit.}
   \item{ols.fit}{LS fitted values.}
   \item{slice.info}{output from sir.slice, used by sir and save.}
   \item{method}{the dimension reduction method used.}
Other returned values repeat quantities from input.  

\code{dr.weights} returns a vector of weights NA substituted for estimated
zero weights.  
}

\references{ The details of these methods are given by R. D. Cook 
(1998).  Regression Graphics.  New York:  Wiley.  Equivalent 
methods are also available in Arc, R. D. Cook and S. Weisberg 
(1999).  Applied Regression Including Computing and Graphics, New 
York:  Wiley, www.stat.umn.edu/arc. } \author{Sanford Weisberg, 
sandy@stat.umn.edu 

For weights, see
R. D. Cook and C. Nachtsheim (1994), Reweighting to achieve
elliptically contoured predictors in regression.  Journal of the American
Statistical Association, 89, 592--599.} %


\seealso{\code{\link{dr.permutation.test}},\code{\link{dr.x}},\code{\link{dr.y}},
\code{\link{dr.direction}},\code{\link{dr.coplot}},\code{\link{dr.weights}}}

\examples{
library(dr)
data(ais)
attach(ais)  # the Australian athletes data
#fit dimension reduction using sir
m1 <- dr(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8)
summary(m1)

# repeat, using save:

m2 <- update(m1,method="save")
summary(m2)

# repeat, using phd:

m3 <- update(m2, method="phdres")
summary(m3)
}

# repeat, using weights:

w1 <- dr.weights(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8)
m4 <- dr(LBM~Wt+Ht+RCC+WCC, method="sir", nslices = 8, weights=w1)
\keyword{regression }%-- one or more ...
