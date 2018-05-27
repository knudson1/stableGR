## used to simulate a multivariate normal using Gibbs



mvn_gibbs <- function(N = 1e4, mu, sigma)
{

  chain <- matrix(0, nrow = N, ncol = p)
  mu1 <- as.matrix(mu[1:(p-1)], ncol = 1)
  mu2 <- mu[p]
  sig11 <- sigma[-p,-p]
  invsig11 <- solve(sig11)
  sig22 <- sigma[p,p]
  invsig22 <- solve(sig22)
  sig12 <- as.matrix(sigma[1:(p-1), p])
  
  # x1 <- as.matrix(rnorm(p-1), ncol = 1)
  # x2 <- rnorm(1)
  x1 <- as.matrix(rep(0, p-1), ncol = 1)
  x2 <- c(0)
  chain[1, ]   <- c(x1, x2)
  for(i in 2:N)
  {
    # Simulating x2 first because will do the same for
    # Linchpin and MWG
    x2mean <- mu2 + t(sig12)%*% invsig11 %*% (x1 - mu1)
    x2sig <- sig22 - t(sig12) %*%invsig11 %*%sig12
    x2 <- rmvnorm(1, mean = x2mean, sigma = x2sig)

    x1mean <- mu1 + sig12 %*% invsig22 %*%t(x2 - mu2)
    x1sig <- sig11 - sig12 %*%invsig22 %*%t(sig12)

    x1 <- t(rmvnorm(1, mean = x1mean, sigma = x1sig))



    chain[i, ]  <- c(x1, x2)
  }

  return(chain)
}
