# Calculates the sqrt(R) = GR diagnostic threshold (arr):
# When the sample diagnostic gets to arr, the chain has converged
# p:dimension of the estimation problem
# m: number of chains
# epsilon: relative precision level
# alpha: significance level

minESS <- mcmcse::minESS

findarr <- function(p, epsilon, m, alpha=.05){
  Tee <- as.numeric(minESS(p, epsilon, alpha = alpha)) #min effective sample size
  del <- sqrt(1 + m/Tee) - 1
  arr <- 1 + del  
  arr
}