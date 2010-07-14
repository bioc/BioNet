\name{hist.bum}
\Rdversion{1.4}
\alias{hist.bum}
\title{
Histogram of the p-value distribution with the fitted bum model
}
\description{
The function plots a histogram of the p-values together with the fitted bum-model.
}
\usage{
\method{hist}{bum}(x, \dots)
}
\arguments{
  \item{x}{
Maximum likelihood estimator object of the beta-uniform mixture fit.}
  \item{\dots}{Other graphic parameters for the plot.}
}
\author{
Daniela Beisser
}
\seealso{\code{\link{fitBumModel}}, \code{\link{hist.bum}}, \code{\link{bumOptim}}}
\examples{
data(pvaluesExample)
pvals <- pvaluesExample[,1]
mle <- fitBumModel(pvals, plot=FALSE)
hist(mle)
}
