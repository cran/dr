\name{fullquad.fit}
\alias{fullquad.fit}
\title{ Internal function }
\description{
This function recovers a linear predictor from a dr object, and then uses
lm to fit a full quadratic regression
}
\usage{
fullquad.fit(object)
}
\arguments{
  \item{object}{ An object of type dr }
}
\value{
  Returns an object of type lm.
}
\author{Sanford Weisberg }

\keyword{internal }% at least one, from doc/KEYWORDS
