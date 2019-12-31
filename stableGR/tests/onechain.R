set.seed(10)
library(mvtnorm)
library(mcmcse)
library(stableGR)

################ 
# Start by making a few chains 
#
# Details on the chain construction
p <- 5
N <- 44^2
tail.ind <- floor(N*.80):N
foo <- matrix(.50, nrow=p, ncol=p)
sigma <- foo^(abs(col(foo)-row(foo)))
mu <- sample(10:20, p)
mu2 <- mu[p]

# Create the chain
mvn.gibbs <- stableGR:::mvn.gibbs
out.gibbs1 <- mvn.gibbs(N = N, mu = mu, sigma = sigma, p = p)
obj <- list(out.gibbs1)

# Convert to MCMC objects
#out1 <- mcmc(out.gibbs1)
#obj <- mcmc.list(out1)

# Calculate psrfs with stable.GR
size <- NULL
method <- "lug"
results1 <- stable.GR(obj, blather = TRUE, size = size, method = method)

##############################################
# Calculate psrfs by for this chain by hand
# Isolate the first component of the  chain
chain1 <- matrix(out.gibbs1[ ,1], ncol = 1)

# Calculate tau^2 for each chain
tausq <- (mcse(chain1, method = method)$se)^2 * N 
tau2 <- asym.var(list(chain1), multivariate = FALSE, method = method, size = size, autoburnin = FALSE)

# Calulate s^2 for each chain
sampvar1 <- var(chain1)
ssquared <- sampvar1 

# Calculate sigma^2 estimate
sigsq <- ((N-1) * ssquared + tausq)/N
othersigsq <- results1$blather$sigmasq
all.equal(as.numeric(sigsq), othersigsq[1])

# Calculate the diagnostic
Rhat <- sigsq / ssquared

that <- sqrt(Rhat)

all.equal(as.numeric(that[1]) , as.numeric(results1$psrf[1]))

