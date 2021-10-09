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
mvn.gibbs <- stableGR:::mvn.gibbs
out.gibbs1 <- mvn.gibbs(N = N, mu = mu, sigma = sigma, p = p)
out.gibbs2 <- mvn.gibbs(N = N, mu = mu, sigma = sigma, p = p)

obj <- list(out.gibbs1, out.gibbs2)

################ 
# Perform unit test using the two chains in obj

outwithfun <- stable.GR(obj, blather = TRUE, size = "sqroot")
withfun <- outwithfun$mpsrf
blather <- outwithfun$blather

# Calculate Tmat for each chain
stacked <- rbind(out.gibbs1, out.gibbs2)
That <- mcse.multi(stacked, method = "lug", size = sqrt(N))$cov
Tmat1 <- mcse.multi(out.gibbs1, method = "lug", size = "sqroot")$cov

all.equal(That, blather$AsymVarMatrix)

#Calc Smat
sloan <- stableGR:::size.and.trim(obj, size = "sqroot")
trimmedchains <- sloan$trimmedchains
stackedchains <- do.call(rbind, trimmedchains)
Smat <- var(stackedchains)
all.equal(Smat, blather$S)


detratio <- det(That)/det(Smat) 
all.equal(detratio, det(solve(blather$S, blather$AsymVarMatrix)))

Nchain <- length(obj)
all.equal(2, Nchain)

top <- (N-1) + ((detratio)^(1/p))
byhand <- sqrt(top/N)

all.equal(byhand, withfun)

################ ################ ################
# Perform unit test using a SINGLE chain (just in case)
onechain <- list(out.gibbs1)
withfun <- stable.GR(onechain, size = "sqroot")$mpsrf

#calculate determinants
Teigen <- eigen(Tmat1)$values
cov1 <- var(out.gibbs1)
Seigen <- eigen(cov1)$values
detT <- (prod(Teigen))
detS <- (prod(Seigen))
detratio <- detT/detS #good

Nchain <- length(onechain)
all.equal(1, Nchain)

rhat <- (N-1)/N + ((detratio)^(1/p))/N
byhand <- sqrt(rhat)
all.equal(byhand, withfun)

