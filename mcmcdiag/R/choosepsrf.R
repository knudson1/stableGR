#' Calculates a Gelman Rubin diagnostic threshold.
#'
#' When the sample diagnostic reaches the psrf threshold calculated in this function, the chain has converged.
#'
#' @param p dimension of the estimation problem
#' @param m number of chains
#' @param delta desired delta value
#' @param epsilon relative precision level
#' @param alpha significance level
#' @return \item{psrf}{The psrf threshold to be reached for convergence.} 
#' @return \item{epsilon}{The epsilon value used to calculate the psrf threshold.}
#' @return \item{delta}{The delta value used to calculate the psrf threshold.}
#' @examples choosepsrf(1, delta = .1, m = 5, alpha = .05)
#' @examples choosepsrf(1, epsilon = .05, m = 5, alpha = .05)


choosepsrf <- function(p, epsilon = .05, delta = NULL, m, alpha=.05){
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





