\name{dr.test}
\alias{dr.test}
\alias{dr.test.default}
\alias{dr.test.phd}
\alias{dr.test.phdq}                        
\alias{dr.test.phdres}                       
\alias{dr.test.phdy}                         
\alias{dr.test.sir}                          
\alias{dr.test2.phdres}    
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Dimension reduction regression tests }
\description{
  This method computes tests of dimension in dimension reduction regression.
  A separate method is needed for each dimension reduction regression method.
  Called by the dr function, and not usually accessed directly by the user.
}
\usage{
dr.test(object, ...) 
dr.test.default(object, nd)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{ dimension reduction regresion object }
  \item{nd}{ Maximum number of dimensions to test }
  \item{...}{arguments passed to dr.test.default or its descendants}
}

\value{
 Returns test statistics and significance levels in a data frame.
}
\references{ R. C. Cook (1998).  Regression Graphics. New York:  Wiley. }
\author{Sanford Weisberg, <sandy@stat.umn.edu> }


\seealso{ \code{\link{dr}} }
\keyword{regression}
