#' Calculates a Gelman Rubin diagnostic threshold.
#'
#' When the sample diagnostic reaches this threshold, the chain has converged.
#'
#' @param p dimension of the estimation problem
#' @param m number of chains
#' @param epsilon relative precision level
#' @param alpha significance level
#' @return The psrf threshold to be reached for convergence.

choosepsrf <- function(p, epsilon, m, alpha=.05){
  Tee <- as.numeric(minESS(p, epsilon, alpha = alpha)) #min effective sample size
  del <- sqrt(1 + m/Tee) - 1
  arr <- 1 + del  
  arr
}

minESS <- mcmcse::minESS





