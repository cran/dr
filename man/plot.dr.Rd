\name{plot.dr}
\alias{plot.dr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Basic plot of a dr object }
\description{
Plots selected direction vectors determined by a dimension reduction regression fit.
By default, the \code{pairs} function is used for plotting, but the user can use any
other graphics command that is appropriate.
}
\usage{
\method{plot}{dr}(x, which = 1:x$numdir, mark.by.y = FALSE, plot.method = pairs, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The name of an object of class dr, a dimension reduction regression object }
  \item{which}{selects the directions to be plotted}
  \item{mark.by.y}{ if TRUE, color points according to the value of the response, otherwise,
  do not color points but include the response as a variable in the plot.}
  \item{plot.method}{the name of a function for the plotting.  The default is \code{pairs}. }
  \item{\dots}{ arguments passed to the plot.method. }
}

\value{
 Returns a graph.
}
\author{Sanford Weisberg, <sandy@stat.umn.edu>.}
 


\examples{
data(ais)
# default fitting method is "sir"
s0 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),data=ais)
plot(s0)
plot(s0,mark.by.y=TRUE)
 }

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }% 
