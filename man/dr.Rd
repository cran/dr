\name{dr}
\alias{dr}
\alias{dr.compute}

\title{Main function for dimension reduction regression}
\description{
 This is the main function in the dr package.  It creates objects of class
 dr to estimate the central (mean) subspace and perform tests concerning
 its dimension.  Several helper functions that require a dr object can then
 be applied to the output from this function.
}
\usage{
dr (formula, data, subset, group=NULL, na.action = na.fail, weights, ...)
    
dr.compute (x, y, weights, group=NULL, method = "sir", chi2approx="bx",...)
 }

\arguments{
 \item{formula}{a two-sided formula like \code{y~x1+x2+x3}, where the left-side
 variable is a vector or a matrix of the response variable(s), and the right-hand side
 variables represent the predictors.  While any legal formula in the Rogers-Wilkinson
 notation can appear, dimension reduction methods generally expect the predictors to be
 numeric, not factors, with no nesting.   Full rank models are
 recommended, although rank deficient models are permitted.
 
 The left-hand side of the formula will generally be a single vector, but it
 can also be a matrix, such as \code{cbind(y1+y2)~x1+x2+x3} if the \code{method}
 is \code{"save"} or \code{"sir"}.  Both of these methods are based on slicing,
 and for the multivariate case slices are determined by slicing on all the
 columns of the left-hand side variables.  
%See \code{\link{dr.slices}} for more information.
 }
 \item{data}{ an optional data frame containing the variables in the model.
           By default the variables are taken from the environment from
           which `dr' is called.}
 \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
 \item{group}{If used, this argument specifies a grouping variable so that 
 dimension reduction is done separately for each distinct level.  This is
 implemented only when \code{method} is one of \code{"sir"}, 
 \code{"save"}, or \code{"ire"}.  This argument must be a one-sided formula.
 For example, \code{~Location} would fit separately for each level of the variable
 \code{Location}.  The formula \code{~A:B} would fit separately for each combination of
 \code{A} and \code{B}, provided that both have been declared factors.}
 \item{weights}{an optional vector of weights to be used where appropriate.  In the
          context of dimension reduction methods, weights are used to obtain
          elliptical symmetry, not constant variance. 
          %See \code{\link{dr.weights}}.
          }
 \item{na.action}{a function which indicates what should happen when the data
          contain `NA's.  The default is `na.fail,' which will stop calculations.
          The option 'na.omit' is also permitted, but it may not work correctly when
          weights are used.}
 \item{x}{The design matrix.  This will be computed from the formula by \code{dr} and then
 passed to \code{dr.compute}, or you can create it yourself.}
 \item{y}{The response vector or matrix}
 \item{method}{This character string specifies the method of fitting.  The options
 include \code{"sir"}, \code{"save"}, \code{"phdy"}, \code{"phdres"} and 
 \code{"ire"}.  Each method may have its own additional arguments, or its
 own defaults; see the details below for more information.}
\item{chi2approx}{Several dr methods compute significance levels using 
     statistics that are asymptotically distributed as a linear combination of
     \eqn{\chi^2(1)} random variables.  This keyword chooses the method for
     computing the chi2approx, either \code{"bx"}, the default for a method
     suggested by Bentler and Xie (2000) or \code{"wood"} for a method proposed
     by Wood (1989).}
 \item{\dots}{For \code{dr}, all additional 
 arguments passed to \code{dr.compute}.  For 
    \code{dr.compute}, additional 
    arguments may be required for particular dimension reduction method.  For
    example, 
    \code{nslices} is the number of slices used by \code{"sir"} and \code{"save"}.
    \code{numdir} is the maximum number of directions to compute, with
    default equal to 4. Other methods may have other defaults.}
}

\details{
The general regression problem studies \eqn{F(y|x)}, the conditional
distribution of a response \eqn{y} given a set of predictors \eqn{x}.  
This function provides methods for estimating the dimension and central
subspace of a general regression problem.  That is, we want to find a 
\eqn{p \times d}{p by d} matrix \eqn{B} of minimal rank \eqn{d} such that 
\deqn{F(y|x)=F(y|B'x)}  
Both the dimension \eqn{d} and the subspace
\eqn{R(B)} are unknown.  These methods make few assumptions.  Many methods
are based on the inverse distribution, \eqn{F(x|y)}.  

For the methods \code{"sir"}, \code{"save"}, \code{"phdy"} and 
\code{"phdres"}, a kernel matrix \eqn{M} is estimated such that the 
column space of \eqn{M} should be close to the central subspace 
\eqn{R(B)}.  The eigenvectors corresponding to the \code{d} largest 
eigenvalues of \eqn{M} provide an estimate of \eqn{R(B)}.

For the method \code{"ire"}, subspaces are estimated by minimizing 
an objective function.

Categorical predictors can be included using the \code{groups} 
argument, with the methods \code{"sir"}, \code{"save"} and 
\code{"ire"}, using the ideas from Chiaromonte, Cook and Li (2002).

The primary output from this method is (1) a set of vectors whose 
span estimates \code{R(B)}; and various tests concerning the 
dimension \code{d}.  

Weights can be used, essentially to specify the relative 
frequency of each case in the data.  Empirical weights that make 
the contours of the weighted sample closer to elliptical can be 
computed using \code{dr.weights}.  
This will usually result in zero weight for some 
cases.  The function will set zero estimated weights to missing.
}

\value{
dr returns an object that inherits from dr (the name of the type is the
value of the \code{method} argument), with attributes:
  \item{x}{The design matrix}
  \item{y}{The response vector}
  \item{weights}{The weights used, normalized to add to n.}
  \item{qr}{QR factorization of x.}
  \item{cases}{Number of cases used.}
  \item{call}{The initial call to \code{dr}.}
  \item{M}{A matrix that depends on the method of computing.  The column space
of M should be close to the central subspace.}
  \item{evalues}{The eigenvalues of M (or squared singular values if M is not
symmetric).}
  \item{evectors}{The eigenvectors of M (or of M'M if M is not square and
symmetric) ordered according to the eigenvalues.}
  \item{chi2approx}{Value of the input argument of this name.}
  \item{numdir}{The maximum number of directions to be found.  The output
value of numdir may be smaller than the input value.}
   \item{slice.info}{output from 'sir.slice', used by sir and save.}
   \item{method}{the dimension reduction method used.} 
   \item{terms}{same as terms attribute in lm or glm.  Needed to make \code{update}
work correctly.}
   \item{A}{If method=\code{"save"}, then \code{A} is a three dimensional array needed to
compute test statistics.}  
}

\references{ 

Bentler, P. M. and Xie, J. (2000), Corrections to test statistics in 
principal Hessian directions.  \emph{Statistics and Probability 
Letters}, 47, 381-389.  Approximate p-values.

Cook, R. D. (1998).  \emph{Regression Graphics}.  New York:  Wiley.  
This book provides the basic results for dimension reduction 
methods, including detailed discussion of the methods \code{"sir"}, 
\code{"phdy"} and \code{"phdres"}.

Cook, R. D. (2004). Testing predictor contributions in sufficient 
dimension reduction. \emph{Annals of Statistics}, 32, 1062-1092.  
Introduced marginal coordinate tests.

Cook, R. D. and Nachtsheim, C. (1994), Reweighting to achieve 
elliptically contoured predictors in regression.  \emph{Journal of 
the American Statistical Association}, 89, 592--599.  Describes the 
weighting scheme used by \code{\link{dr.weights}}.

Cook, R. D. and Ni, L. (2004). Sufficient dimension reduction via 
inverse regression:  A minimum discrrepancy approach, \emph{Journal 
of the American Statistical Association}, 100, 410-428. The 
\code{"ire"} is described in this paper.

Cook, R. D. and Weisberg, S. (1999).  \emph{Applied Regression 
Including Computing and Graphics}, New York:  Wiley, 
\url{http://www.stat.umn.edu/arc}.  The program \code{arc} described 
in this book also computes most of the dimension reduction methods 
described here.

Chiaromonte, F., Cook, R. D. and Li, B. (2002). Sufficient dimension 
reduction in regressions with categorical predictors. Ann. Statist. 
30 475-497.  Introduced grouping, or conditioning on factors.
 
Shao, Y., Cook, R. D. and Weisberg (2007).  Marginal tests with 
sliced average variance estimation.  \emph{Biometrika}.  Describes 
the tests used for \code{"save"}. 

Wen, X. and Cook, R. D. (2007).  Optimal Sufficient Dimension 
Reduction in Regressions with Categorical Predictors, \emph{Journal 
of Statistical Inference and Planning}.   This paper extends the 
\code{"ire"} method to grouping.  

Wood, A. T. A. (1989) An \eqn{F} approximation to the distribution 
of a linear combination of chi-squared variables. 
\emph{Communications in Statistics: Simulation and Computation}, 18, 
1439-1456.  Approximations for p-values. } 

\author{Sanford Weisberg, <sandy@stat.umn.edu>.} %


%\seealso{\code{\link{dr.permutation.test}},\code{\link{dr.coordinate.test}},
%\code{\link{dr.direction}},\code{\link{dr.coplot}},\code{\link{dr.weights}},
%\code{\link{dr.pvalue}}}

\examples{
data(ais)
# default fitting method is "sir"
s0 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),data=ais) 
# Refit, using a different function for slicing to agree with arc.
summary(s1 <- update(s0,slice.function=dr.slices.arc))
# Refit again, using save, with 10 slices; the default is max(8,ncol+3)
summary(s2<-update(s1,nslices=10,method="save"))
# Refit, using phdres.  Tests are different for phd, and not
# Fit using phdres; output is similar for phdy, but tests are not justifiable. 
summary(s3<- update(s1,method="phdres"))
# fit using ire:
summary(s4 <- update(s1,method="ire"))
# fit using Sex as a grouping variable.  
s5 <- update(s4,group=~Sex)
}

\keyword{regression }%-- one or more ...
