set.seed(10)
library(mvtnorm)
library(mcmcse)
library(coda)
library(mcmcdiag)

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
mvn_gibbs <- mcmcdiag:::mvn_gibbs
out.gibbs1 <- mvn_gibbs(N = N, mu = mu, sigma = sigma, p = p)
out.gibbs2 <- mvn_gibbs(N = N, mu = mu, sigma = sigma, p = p)

# Convert to MCMC objects
out1 <- mcmc(out.gibbs1)
out2 <- mcmc(out.gibbs2)
obj <- mcmc.list(out1, out2)

################ 
# Perform unit test using the two chains in obj

outwithfun <- gelman.bm(obj, blather = TRUE)
withfun <- outwithfun$mpsrf
blather <- outwithfun$blather

# Calculate Tmat for each chain
Tmat1 <- mcse.multi(out1, method = "lug")$cov
Tmat2 <- mcse.multi(out2, method = "lug")$cov
That <- .5*(Tmat1 + Tmat2) 
all.equal(That, blather$Tee)

#Calc Smat
cov1 <- var(out1)
cov2 <- var(out2)
Smat <- .5*(cov1 + cov2) 
all.equal(Smat, blather$S)

#calculate determinants
#Teigen <- eigen(That)$values
#Seigen <- eigen(Smat)$values
#detT <- (prod(Teigen))
#detS <- (prod(Seigen))
detratio <- det(That)/det(Smat) 
all.equal(detratio, det(solve(blather$S, blather$Tee)))

Nchain <- nchain(obj)
all.equal(2, Nchain)

top <- (N-1) + ((detratio)^(1/p))
byhand <- sqrt(top/N)

all.equal(byhand, withfun)

################ ################ ################
# Perform unit test using a SINGLE chain (just in case)
onechain <- mcmc.list(out1)
withfun <- gelman.bm(onechain)$mpsrf

#calculate determinants
Teigen <- eigen(Tmat1)$values
Seigen <- eigen(cov1)$values
detT <- (prod(Teigen))
detS <- (prod(Seigen))
detratio <- detT/detS #good

Nchain <- nchain(onechain)
all.equal(1, Nchain)

rhat <- (N-1)/N + ((detratio)^(1/p))/N
byhand <- sqrt(rhat)
all.equal(byhand, withfun)

