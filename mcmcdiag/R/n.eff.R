#' Effective sample size 
#' 
#' For an estimator, effective sample size is the number of independent samples with the same standard error as a correlated sample. This function calculates effective sample size for a set of Markov chains using lugsail variance estimators. This also determines whether the Markov chains have converged. If they have not, this function approximates the number of samples needed.
#'
#' @param x an \code{mcmc.list} object with more than one chain,  with starting values that are overdispersed with respect to the posterior distribution.
#' @param ... arguments to be passed to \code{gr.diag}.
#'
#' @return \item{n.eff}{a scalar point estimate of the effective sample size.}
#' @return \item{converged}{a logical indicating whether the sample has converged.}
#' @return \item{n.target}{NULL (if the chain has converged) or a scalar estimate of the chain length required for convergence, assuming the number of chains is unchanged.  }
#' @return \item{n.more}{NULL (if the chain has converged) or a scalar estimate of number of additional samples required for convergence (per chain), assuming the number of chains is unchanged.  }
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
  
  # prepare and do comparison to our goal
  currentESS <- out$n.eff
  p <- nvar(x)
  m <- nchain(x)
  targ <- target.psrf(p=p, m=m)
  targetESS <- RtoESS(targ$psrf, m)
  converged <- FALSE
  if(currentESS >= targetESS){ converged <- TRUE}
  
  # if the sample hasn't converged, approximate the sample size that would result in convergence
  ntarget <- nmore <- NULL
  if(converged == FALSE){
    ncurrent <- niter(x)
    ntarget <- ceiling(ncurrent*targetESS/currentESS)
    nmore <- ceiling(ntarget - ncurrent)
  }  
  
  list(n.eff = currentESS, converged = converged, n.target = ntarget, n.more = nmore)
}

RtoESS <- function(arrr,  m){
  m/(arrr^2 - 1)
}


