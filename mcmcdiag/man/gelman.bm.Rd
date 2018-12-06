\name{gelman.bm}
\alias{gelman.bm}

\title{A Gelman-Rubin diagnostic with batch means}

\description{This function uses batch means estimators to calculate a convergence diagnostic for Markov chain Monte Carlo in the spirit of Gelman-Rubin.

 A univariate `potential scale reduction factor' (PSRF) is calculated for each variable in
\code{x}. For multivariate chains, a multivariate PSRF is calculated to take into account the interdependence of the chain's components.  The psrfs decreases to 1 as the chain length increases. When the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has converged to to the target distribution  onvergence.
}

\usage{
gelman.bm(x, mapping = "determinant",  multivariate = TRUE, 
	method = "lug", blather = FALSE)
}

\arguments{
\item{x}{an \code{mcmc.list} object with more than one chain,
and with starting values that are overdispersed with respect to the
posterior distribution.}
\item{mapping}{the function used to map the covariance matrix to a scalar. This is one of 
\dQuote{\code{determinant}} (determinant of the covariance matrix, 
  the default) or \dQuote{\code{maxeigen}} (the largest eigenvalue of the covariance matrix).}

\item{multivariate}{a logical flag indicating whether the multivariate
potential scale reduction factor should be calculated for multivariate
chains}
 \item{method}{the method used to compute the standard
  error of the chains. This is one of \dQuote{\code{bm}} (batch means,
  the default), \dQuote{\code{obm}} (overlapping batch
  means), \dQuote{\code{tukey}} (spectral variance method
  with a Tukey-Hanning window), or \dQuote{\code{bartlett}}
  (spectral variance method with a Bartlett window).}

\item{blather}{a logical flag indicating whether to include additional output.}

}



\value{

  An object of class \code{gelman.diag}. This is a list with the
  following elements:
  \item{psrf}{A vector containing the point estimates of the potential
  scale reduction factor.}
    \item{mpsrf}{A scalar point estimate of the multivariate potential
  scale reduction factor.}
      \item{means}{A vector containing the sample means based on the chains provided.}
	\item{blather}{Either \code{FALSE} or a list containing intermediate calculations.}
}


\section{Theory}{ The univariate  potential scale reduction factor (PSRF) was first Gelman and Rubin (1992)

The univariate PSRF is the square root of the ratio of two estimators for the target variance. The PSRF decreases to 1 as the chain length increases; when the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has converged to to the target distribution.}


\references{
  Flegal, J. M. and Jones, G. L. (2010) Batch means and
  spectral variance estimators in Markov chain Monte Carlo.
  \emph{The Annals of Statistics}, \bold{38}, 1034--1070.

Gelman, A and Rubin, DB (1992) Inference from iterative simulation
using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511.

Brooks, SP. and Gelman, A. (1998) General methods for monitoring
convergence of iterative simulations. \emph{Journal of Computational and
Graphical Statistics}, \bold{7}, 434-455.
}

\author{
Christina Knudson
}

 
\keyword{  htest }

