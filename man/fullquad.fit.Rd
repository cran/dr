\name{fullquad.fit}
\alias{fullquad.fit}
\alias{fullquad.residuals}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Fit a full second-order model}
\description{
  This function fits an ordinary LS regression of a response on a full
  second order model, or returns the residuals from such a fit. Used in the
  phdq method in dr.
}
\usage{
fullquad.fit(x, y, w, ...)
fullquad.residuals(...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{ an n by p matrix of predictors}
  \item{y}{ an n vector of the response}
  \item{w}{ an n vector of positive weights or NULL}
  \item{\dots}{ arguments passed to lm}
}
\value{
  fullquad.fit returns an object of type lm.  fullquad.residuals returns a
  vector of n residuals.
}
\author{ Sanford Weisberg <sandy@stat.umn.edu>}



\keyword{ regression }% at least one, from doc/KEYWORDS
