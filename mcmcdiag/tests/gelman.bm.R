set.seed(10)
library(mvtnorm)
library(mcmcse)

################ 
#Start by making a few chains to work with
#Details on the chain construction
mvn_gibbs <- mcmcdiag:::mvn_gibbs
p <- 5
N <- 5e2
tail.ind <- floor(N*.80):N
foo <- matrix(.50, nrow=p, ncol=p)
sigma <- foo^(abs(col(foo)-row(foo)))
mu <- sample(10:20, p)
mu2 <- mu[p]
#Create the chains
out.gibbs1 <- mvn_gibbs(N = N, mu = mu, sigma = sigma)
out.gibbs2 <- mvn_gibbs(N = N, mu = mu, sigma = sigma)
#Convert to MCMC objects
out1 <- mcmc(out.gibbs1)
out2 <- mcmc(out.gibbs2)
obj <- mcmc.list(out1, out2)

################ 
gelman.bm(obj)



outmcse <- mcse.mat(out.gibbs1)
outmcse

# calculating taus
this <- sapply(x, function(x1) {(mcse.mat(x1)[ ,2])^2 *Niter})
that <- cbind((mcse.mat(out.gibbs1)[ ,2])^2 * Niter, (mcse.mat(out.gibbs2)[ ,2])^2 * Niter)
all.equal(this, that) 
# ncol = number of chains.
# nrow = number of variables

 taus <- matrix(sapply(x, apply, 3, mcse.mat, simplify = TRUE), 
        nrow = Nvar, ncol = Nchain) #calculates column means for each chain



