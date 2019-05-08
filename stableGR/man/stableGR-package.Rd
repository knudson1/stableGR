\name{stableGR-package}
\alias{stableGR-package}
\alias{stableGR}
\docType{package}
\title{
\packageTitle{stableGR}
}
\description{
% Uses stable variance estimates to calculate the Gelman-Rubin diagnostic, in addition to smart cut-offs for the diagnostic.
\packageDescription{stableGR}
}
\details{

\packageDESCRIPTION{stableGR}
% \packageIndices{stableGR}
The differentiating feature of this package is the connection between effective sample size and the potential scale reduction factor (PSRF), compatible with multiple chains. The two main functions in the package are \code{stable.GR}, \code{n.eff}, and \code{target.psrf}. \code{stable.GR} will return the univariate and multivariate PSRF, along with an estimate of the effective sample size. \code{n.eff} will return an estimate of the effective sample size, along with an estimate of the target effective sample size that should be reached. \code{target.psrf} returns the desired cutoff for \code{stable.GR}, so that simulation of MCMC should stop only when the \code{target.psrf} values has been obtained.
}
\author{
\packageAuthor{stableGR}

Maintainer: \packageMaintainer{stableGR}
}
\references{ 
Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384. 

Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.

Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. 

Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. 

Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
}


% \keyword{ htest }
\examples{
library(stableGR)
set.seed(100)
p <- 2
n <- 1000

sig.mat = matrix(c(1, .8, .8, 1), ncol = 2, nrow = 2)

# Making 3 chains
\dontrun{chain1 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
chain2 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
chain3 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)

# find GR diagnostic using all three chains
x <- list(chain1, chain2, chain3)
stable.GR(x) 
}

}
