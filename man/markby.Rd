\name{markby}
\alias{markby}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Produce a vector of colors/symbols to mark points}
\description{
This function creates marking of points by color or symbol for use in graphs.
}
\usage{
markby(z, use="color", values=NULL, color.fn=rainbow, na.action="na.use")
}
\arguments{
 \item{z}{a variable with a few distict values that will define the groups for
marking.}
 \item{use}{if equal to \code{"color"}, will returna list of colors.  If
anything else, it will return a list of symbols for marking.}
 \item{values}{a list with as many values as unique values of z that determine
the colors or symbols.  If this is not set, then the function rainbow is used
for colors and 1:length(unique(z)) is used for symbols.}
 \item{color.fn} {This is the name of a function that when called with a
 single argument of length k returns a list of k acceptable color names.
 The function rainbow is a standard function in \R; see the help for this
 function for other choices.}
 \item{na.action}{By default, missing values (NAs) are treated as another
 category in the color marking.  Setting this argument to anything else will
 return NA for the missing category names.  This will cause an error in \R\
 graphics routines, but see below for usage of this options.}
}
\value{
Returns length(z) values that specify the color or symbol for each point.
}
\references{This function is to help users familiar with Arc, as discussed in
R. D. Cook and S. Weisberg (1999).  \emph{Applied Regression Including
Computing and Graphics}, New York:  Wiley.}
\author{Sanford Weisberg, <sandy@stat.umn.edu>}


\examples{
  x <- rnorm(100)
  y <- rnorm(100)
  z <- cut(rnorm(100),3)
# Scatterplot, mark using color with groups determined by Status
  plot(x,y,col=markby(z))
# Scatterplot, mark using symbols with groups determined by Status
  plot(x,y,pch=markby(z,use="symbols"))
# handing of missing values:
  z[1:10] <- NA   # set to missing
  marks <- markby(z,na.action="omit")
  sel <- !is.na(marks)
  plot(x[sel],y[sel],col=marks[sel])
}
\keyword{ regression }%-- one or more ...
