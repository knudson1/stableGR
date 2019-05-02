#' Two block Gibbs sampler for a multivariate normal distribution 
#' 
#' This function generates a Markov chain sample from a multivariate normal distribution using a two-block Gibbs sampler. The function is used mainly for implementation in the examples of this package.
#' 
#' @param N number of Markov chain samples desired
#' @param p dimension of the multivariate normal target distribution
#' @param mu mean vector of the multivariate normal distribution
#' @param sigma covariance matrix of the multivariate normal distribution
#'
#' @return  N by p matrix of samples from the multivariate normal target distribution
#'
#' @export
#'

mvn.gibbs <- function(N = 1e4,p, mu, sigma)
{

  chain <- matrix(0, nrow = N, ncol = p)
  mu1 <- as.matrix(mu[1:(p-1)], ncol = 1)  # first block
  mu2 <- mu[p]    #second block
  sig11 <- sigma[-p,-p]
  invsig11 <- solve(sig11)
  sig22 <- sigma[p,p]
  invsig22 <- solve(sig22)
  sig12 <- as.matrix(sigma[1:(p-1), p])
  

  x1 <- as.matrix(rep(0, p-1), ncol = 1) # starting value is zero
  x2 <- c(0)  #zero starting value
  chain[1, ]   <- c(x1, x2)
  for(i in 2:N)
  {

    x2mean <- mu2 + t(sig12)%*% invsig11 %*% (x1 - mu1)
    x2sig <- sig22 - t(sig12) %*%invsig11 %*%sig12
    x2 <- rmvnorm(1, mean = x2mean, sigma = x2sig)  #sample x2

    x1mean <- mu1 + sig12 %*% invsig22 %*%t(x2 - mu2)
    x1sig <- sig11 - sig12 %*%invsig22 %*%t(sig12)

    x1 <- t(rmvnorm(1, mean = x1mean, sigma = x1sig))  #sample x1
    chain[i, ]  <- c(x1, x2)
  }

  return(chain)
}
