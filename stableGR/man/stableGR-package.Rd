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
Practitioners of Bayesian statistics often use Markov chain Monte Carlo (MCMC) samplers to sample from a posterior distribution. This package determines whether the MCMC sample is large enough   to yield reliable estimates of the target distribution. In particular, this calculates a Gelman-Rubin convergence diagnostic using stable and consistent estimators of Monte Carlo variance. Additionally, this uses the connection between an MCMC sample's effective sample size and the Gelman-Rubin diagnostic to produce a threshold for terminating MCMC simulation. Finally, this informs the user whether enough samples have been collected  and (if necessary) estimates the number of samples needed for a desired level of accuracy. The theory underlying these methods can be found in "Revisiting the Gelman-Rubin Diagnostic" by Vats and  Knudson (2018) <arXiv:1812:09384>. 


This package is unique in a few ways. First, it uses stable variance estimators to calculate a stabilized Gelman-Rubin statistic. Second, it leverages the connection between effective sample size and the potential scale reduction factor (PSRF). Third, this diagnostic can be used whether MCMC samples were created from a single chain or multiple chains. 

The main functions in the package are \code{stable.GR}, \code{n.eff}, and \code{target.psrf}. \code{stable.GR}  returns the univariate PSRF, the multivariate PSRF, and the estimated effective sample size. \code{n.eff}  returns informs the user whether sufficient MCMC samples have been collected; if not, \code{n.eff} also returns the estimated target sample size \code{target.psrf} creates a termination threshold for \code{stable.GR}; MCMC sampling can terminate when the MCMC samples' psrf is smaller than the value returned by \code{target.psrf}.
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
n <- 100 # For real problems, use a MUCH larger n.

sig.mat = matrix(c(1, .8, .8, 1), ncol = 2, nrow = 2)

# Making 3 chains
chain1 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
chain2 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
chain3 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)

# find GR diagnostic using all three chains
x <- list(chain1, chain2, chain3)
stable.GR(x) 


}
