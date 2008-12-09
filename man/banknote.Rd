\name{banknote}
\alias{banknote}
\non_function{}
\title{Swiss banknote data}
\description{
Six measurements made on 100 genuine
Swiss banknotes and 100 counterfeit ones.
}
\format{This data frame contains the following columns:
 \describe{
    \item{Length}{
 Length of bill, mm
    }
    \item{Left}{
 Width of left edge, mm
    }
    \item{Right}{
 Width of right edge, mm
    }
    \item{Bottom}{
 Bottom margin width, mm
    }
    \item{Top}{
 Top margin width, mm
    }
    \item{Diagonal}{
 Length of image diagonal, mm
    }
    \item{Y}{
 0 = genuine, 1 = counterfeit
    }
 }
}

\source{
Flury, B. and Riedwyl, H. (1988). \emph{Multivariate Statistics: A practical approach.} London: Chapman & Hall.}
\references{Weisberg, S. (2005).  \emph{Applied Linear Regression}, 3rd edition. New York: Wiley, Problem 12.5.}
\examples{
data(banknote)
}
\keyword{datasets}
