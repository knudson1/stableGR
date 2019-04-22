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

obj <- list(out.gibbs1, out.gibbs2)

# Convert to MCMC objects
#out1 <- mcmc(out.gibbs1)
#out2 <- mcmc(out.gibbs2)
#obj <- mcmc.list(out1, out2)

################ ################ 
# Perform unit test using the two chains in obj

outwithfun <- gr.diag(obj, blather = TRUE)
mpsrf <- outwithfun$mpsrf
blather <- outwithfun$blather
m <- length(obj)

denom <- mpsrf^2 - ((N-1)/N)
byhand <- m/denom

all.equal(byhand, outwithfun$n.eff)
all.equal(byhand, as.numeric(n.eff(obj)[1]))
 
################ ################ 
# make sure arguments can be passed through n.eff to gr.diag correctly
################ 
outA <- gr.diag(obj, autoburnin = TRUE)
outB <- n.eff(obj, autoburnin = TRUE)
all.equal(outA$n.eff, outB$n.eff)

################ 
outA <- gr.diag(obj, method = "bm")
outB <- n.eff(obj, method = "bm")
all.equal(outA$n.eff, outB$n.eff)

################ ################ 
# make sure this works for univariate
babyout1 <- out1[,1]
babyout2 <- out2[,1]
babyobj <- mcmc.list(babyout1, babyout2)
temp <- gr.diag(babyobj)$psrf
denom <- temp^2 - ((N-1)/N)
m <- length(babyobj)
byhand <- m/denom
melon <- n.eff(babyobj)
all.equal(melon$n.eff, byhand)
melon$converged

################ ################ 
set.seed(1234)
banana <- mcmc(rnorm(1000))
banana <- mcmc.list(banana)
n.eff(banana)$converged == FALSE
