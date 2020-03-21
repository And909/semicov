\name{semCov}
\alias{semCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Approximated Semicovariance Matrix
}
\description{
This function computes an approximated semicovariance matrix as in Estrada (2008).
}
\usage{
semCov(mydata,bench)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata}{Matrix of returns. Each row is a period and each column is an asset.}
  \item{bench}{Benchmark for semivariance. Defaults to 0.}
}
\details{
See Rigamonti (2020) for more details.
}
\value{
A matrix that approximates the semicovariance matrix.
}
\references{
Estrada, J. (2008). Mean-semivariance optimization: A heuristic approach. Journal of Applied Finance, 18(1):57-72
Rigamonti, A. (2020). Numerical solution for the minimum semivariance portfolio optimization problem in R: the semicov package
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.

%%## The function is currently defined as
%%function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
