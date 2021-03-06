% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvn.gibbs.R
\name{mvn.gibbs}
\alias{mvn.gibbs}
\title{Two block Gibbs sampler for a multivariate normal distribution}
\usage{
mvn.gibbs(N = 10000, p, mu, sigma)
}
\arguments{
\item{N}{number of Markov chain samples desired}

\item{p}{dimension of the multivariate normal target distribution}

\item{mu}{mean vector of the multivariate normal distribution}

\item{sigma}{covariance matrix of the multivariate normal distribution}
}
\value{
N by p matrix of samples from the multivariate normal target distribution
}
\description{
This function generates a Markov chain sample from a multivariate normal distribution using a two-block Gibbs sampler. The function is used mainly for implementation in the examples of this package.
}
