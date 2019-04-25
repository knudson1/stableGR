set.seed(10)
library(mvtnorm)
library(mcmcse)
library(stableGR)

################ 
# Start by making a few chains to work with

# Details on the chain construction
p <- 5
N <- 10000
tail.ind <- floor(N*.80):N
foo <- matrix(.50, nrow=p, ncol=p)
sigma <- foo^(abs(col(foo)-row(foo)))
mu <- sample(10:20, p)
mu2 <- mu[p]

# Create the chains
mvn_gibbs <- stableGR:::mvn_gibbs
out.gibbs1 <- mvn_gibbs(N = N, mu = mu, sigma = sigma, p = p)
out.gibbs2 <- mvn_gibbs(N = N, mu = mu, sigma = sigma, p = p)

obj <- list(out.gibbs1, out.gibbs2)

# Convert to MCMC objects
#out1 <- mcmc(out.gibbs1)
#out2 <- mcmc(out.gibbs2)
#obj <- mcmc.list(out1, out2)

################ 
# Perform unit test using the two chains in obj

withfun <- stable.GR(obj, mapping = "maxeigen")$mpsrf


# Calculate Tmat for each chain
Tmat1 <- mcse.multi(out.gibbs1, method = "lug")$cov
Tmat2 <- mcse.multi(out.gibbs2, method = "lug")$cov
That <- .5*(Tmat1 + Tmat2) #good

#Calc Smat
cov1 <- var(out.gibbs1)
cov2 <- var(out.gibbs2)
Smat <- .5*(cov1 + cov2) #good

Sinv <- qr.solve(Smat)
thingy <- Sinv %*% That
maxeigen <- max(eigen(thingy)$values)

Nchain <- length(obj)
all.equal(2, Nchain)

rhat <- (N-1)/N +  maxeigen/N
byhand <- sqrt(rhat)
all.equal(byhand, withfun)
