\name{cosangle}
\alias{cosangle}
\alias{cosangle1}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{ Compute the cosine of the angle between a vector and a subspace }
\description{
\code{cosangle1} returns the cosine of the
angle between a vector vecs  and the subspace spanned
by the columns of the matrix mat.  For cosangle, vecs can be a matrix, in
which case cosangle1 is called for each column in vecs.
}
\usage{
cosangle(mat, vecs)
cosangle1(mat, vec)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mat}{ A matrix the provides the basis for a subspace}
  \item{vecs}{ A vector with the same number of rows as mat, or a matrix
  with the same number of rows as mat}
  \item{vec}{ A vector with the same number of rows as mat}
}
\details{
  \code{cosangle1} computes the cosine of the angle between vec and the orthogonal projection
  of vec onto the column space of mat.  \code{cosangle} repeats this computation
  for each column of vecs.
}
\value{
  cosangle1 returns a single value and cosangle returns a vector.
}
\author{Sanford Weisberg <sandy@stat.umn.edu> }
\keyword{regression}
