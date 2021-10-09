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


################ ################ 
# Perform unit test using the two chains in obj

outwithfun <- stable.GR(obj, blather = TRUE)
mpsrf <- outwithfun$mpsrf
blather <- outwithfun$blather
m <- length(obj)
Nneeded <- blather$Nneeded

denom <- mpsrf^2 - ((N-1)/N)
denom <- mpsrf^2 - ((Nneeded-1)/Nneeded)
byhand <- m/denom

all.equal(byhand, outwithfun$n.eff)
all.equal(byhand, as.numeric(n.eff(obj)[1]))
 
################ ################ 
# make sure arguments can be passed through n.eff to stable.GR correctly
################ 
outA <- stable.GR(obj, autoburnin = FALSE)
outB <- n.eff(obj, autoburnin = FALSE)
all.equal(outA$n.eff, outB$n.eff)

################ 
outA <- stable.GR(obj, method = "bm")
outB <- n.eff(obj, method = "bm")
all.equal(outA$n.eff, outB$n.eff)

################ ################ 
# make sure this works for univariate
babyout1 <- matrix(out.gibbs1[,1], ncol = 1)
babyout2 <- matrix(out.gibbs2[,1], ncol = 1)
babyobj <- list(babyout1, babyout2)
suppressWarnings(confrontyourracistuncle <- stable.GR(babyobj)$psrf)
denom <- confrontyourracistuncle^2 - ((N-1)/N)
m <- length(babyobj)
byhand <- m/denom
suppressWarnings(melon <- n.eff(babyobj))
all.equal(melon$n.eff, byhand)
melon$converged

################ ################ 
set.seed(1234)
banana <- list(matrix(rnorm(1000), ncol = 1))
suppressWarnings(n.eff(banana)$converged) == FALSE

