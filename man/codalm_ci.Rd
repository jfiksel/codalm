% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codalm_em.R
\name{codalm_ci}
\alias{codalm_ci}
\title{Bootstrap Confidence Intervals Linear Regression for Compositional Outcomes and Predictors}
\usage{
codalm_ci(
  y,
  x,
  accelerate = TRUE,
  nboot = 500,
  conf = 0.95,
  parallel = FALSE,
  ncpus = NULL,
  strategy = NULL,
  init.seed = 123
)
}
\arguments{
\item{y}{A matrix of compositional outcomes. Each row is an observation, and must sum to 1.
If any rows do not sum to 1, they will be renormalized}

\item{x}{A matrix of compositional predictors. Each row is an observation, and must sum to 1.
If any rows do not sum to 1, they will be renormalized}

\item{accelerate}{A logical variable, indicating whether or not to use the
Squarem algorithm for acceleration of the EM algorithm. Default is TRUE}

\item{nboot}{The number of bootstrap repetitions to use. Default is 500}

\item{conf}{A scalar between 0 and 1 containing the confidence level of the required intervals.
Default is .95.}

\item{parallel}{A logical variable, indicating whether or not to use a parallel
operation for computing the permutation statistics}

\item{ncpus}{Optional argument. When provided, is an integer giving the number
of clusters to be used in parallelization. Defaults to the number of cores, minus 1.}

\item{strategy}{Optional argument. When provided, this will be the evaluation function
(or name of it) to use for parallel computation (if parallel = TRUE). Otherwise,
if parallel = TRUE, then this will default to multisession. See \code{\link[future]{plan}}.}

\item{init.seed}{The initial seed for the permutations. Default is 123.}
}
\value{
A list, with \code{ci_L} and \code{ci_U}, giving the lower and upper bounds
of each element of the B matrix
}
\description{
Implements percentile based bootstrapping to estimate the confidence intervals
for the regression coefficients when doing linear regression for compositional outcomes
and predictors
}
\examples{
\donttest{
require(ggtern)
data("WhiteCells", package = 'ggtern')
image <- subset(WhiteCells, Experiment == "ImageAnalysis")
image_mat <- as.matrix(image[,c("G", "L", "M")])
microscopic <- subset(WhiteCells, Experiment == "MicroscopicInspection")
microscopic_mat <- as.matrix(microscopic[,c("G", "L", "M")])
x <- image_mat  / rowSums(image_mat)
y <- microscopic_mat / rowSums(microscopic_mat)
codalm_ci(y, x, nboot = 50, conf = .95)
}
}
