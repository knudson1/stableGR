set.seed(10)
library(mvtnorm)
library(mcmcse)
library(coda)

################ 
# Start by making a few chains to work with
# Details on the chain construction
mvn_gibbs <- mcmcdiag:::mvn_gibbs
p <- 5
N <- 100
tail.ind <- floor(N*.80):N
foo <- matrix(.50, nrow=p, ncol=p)
sigma <- foo^(abs(col(foo)-row(foo)))
mu <- sample(10:20, p)
mu2 <- mu[p]
# Create the chains
out.gibbs1 <- mvn_gibbs(N = N, mu = mu, sigma = sigma)
out.gibbs2 <- mvn_gibbs(N = N, mu = mu, sigma = sigma)
# Convert to MCMC objects
out1 <- mcmc(out.gibbs1)
out2 <- mcmc(out.gibbs2)
obj <- mcmc.list(out1, out2)

################ 
# Perform unit test using the two chains in obj
this <- gelman.bm(obj)

#all.equal(this$psrf, that)



