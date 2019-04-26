\name{stableGR-package}
\alias{stableGR-package}
\alias{stableGR}
\docType{package}
\title{
\packageTitle{stableGR}
}
\description{Uses batch means to calculate a convergence diagnostic for Markov chain Monte Carlo in the spirit of Gelman-Rubin
\packageDescription{stableGR}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{stableGR}
\packageIndices{stableGR}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
\packageAuthor{stableGR}

Maintainer: \packageMaintainer{stableGR}
}
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


\keyword{ htest }
\examples{
library(stableGR)
set.seed(100)

p <- 5
n <- 10000

# Making 3 chains
chain1 <- matrix(rnorm(p*n), ncol = p, nrow = n)
chain2 <- matrix(rnorm(p*n), ncol = p, nrow = n)
chain3 <- matrix(rnorm(p*n), ncol = p, nrow = n)

# find GR diagnostic using all three chains
x <- list(chain1, chain2, chain3)
stable.GR(x) 

# find GR diagnostic using a single chain
stable.GR(chain1)


}
