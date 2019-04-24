#' Calculates a Gelman Rubin diagnostic threshold using lugsail variance estimators.
#'
#' When the sample diagnostic reaches the psrf threshold calculated in this function, the chain has converged.
#'
#' @param p dimension of the estimation problem.
#' @param epsilon relative precision level.
#' @param delta desired delta value.
#' @param m number of chains.
#' @param alpha significance level.
#' 
#' @return \item{psrf}{The psrf threshold to be reached for convergence.} 
#' @return \item{epsilon}{The epsilon value used to calculate the psrf threshold.}
#' @return \item{delta}{The delta value used to calculate the psrf threshold.}
#' @examples target.psrf(1, delta = .1, m = 5, alpha = .05)
#' @examples target.psrf(1, epsilon = .05, m = 5, alpha = .05)
#' 
#' @section References:
#' Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384 
#' 
#' Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.
#' 
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. \cr
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. \cr
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @export
#'

target.psrf <- function(p, epsilon = .05, delta = NULL, m, alpha=.05){
  if(is.null(delta)){
    Tee <- as.numeric(minESS(p, epsilon, alpha = alpha)) #min effective sample size
  }
  
  if(is.null(delta) == FALSE){
    Tee <- M <- m/((1+delta)^2 - 1)
    epsilon <- as.numeric(minESS(p, eps = epsilon, ess = M))
  }
  
  del <- sqrt(1 + m/Tee) - 1
  arr <- 1 + del  
  list(psrf = arr, epsilon = epsilon, delta = del)
}

minESS <- mcmcse::minESS





