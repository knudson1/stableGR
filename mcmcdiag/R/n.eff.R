#' Effective sample size 
#' 
#' For an estimator, effective sample size is the number of independent samples with the same standard error as a correlated sample. This function calculates effective sample size for a set of Markov chains using lugsail variance estimators. This also determines whether the Markov chains have converged. If they have not, this function approximates the number of samples needed.
#'
#' @param x an \code{mcmc.list} object with more than one chain,  with starting values that are overdispersed with respect to the posterior distribution.
#' @param ... arguments to be passed to \code{gr.diag}.
#'
#' @return \item{n.eff}{A scalar point estimate of the effective sample size.}
#'
#'
#' @section References:
#' Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384. 
#' 
#' Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.
#' 
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. \cr
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. \cr
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @export
#'


n.eff <- function(x, ...){ 
  out <- gr.diag(x, ...)
  list(n.eff = out$n.eff)
}

gr.diag <- mcmcdiag:::gr.diag
