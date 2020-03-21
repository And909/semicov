\name{numMinSvar}
\alias{numMinSvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum Semivariance Weights
}
\description{
This function numerically computes the weights for the minimum semivariance portfolio with up to five assets.
}
\usage{
numMinSvar(mydata,ngrid,bench,short)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mydata}{Matrix of returns. Each row is a period and each column is an asset.}
  \item{ngrid}{Number of values in the grid of weights. Defaults to 201.}
  \item{bench}{Benchmark for semivariance. Defaults to 0.}
  \item{short}{Allowing short selling. Defaults to FALSE (no short positions allowed).}
}
\details{
Short selling allowed with up to 3 assets. Weights range with short selling spans from -0.5 to 1.5.
Maximum ngrid allowed: 201 with short selling and 3 assets; 201 with 4 assets; 101 with 5 assets.
ngrid can be raised above the default level of 201 when using 2 assets or 3 assets without short selling. 
Maximum number of periods allowed in the estimation window with 5 assets is 120.
See Rigamonti (2020) for more details.
}
\value{
A numeric vector with the minimum semivariance portfolio weights.
}
\references{
Rigamonti, A. (2020). Numerical solution for the minimum semivariance portfolio optimization problem in R: the semicov package
}
\author{
Andrea Rigamonti
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
