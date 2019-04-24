set.seed(10)
library(mvtnorm)
library(mcmcse)
library(stableGR)

################ 
# Start by making a few chains 
#
# Details on the chain construction
p <- 5
N <- 2000
tail.ind <- floor(N*.80):N
foo <- matrix(.50, nrow=p, ncol=p)
sigma <- foo^(abs(col(foo)-row(foo)))
mu <- sample(10:20, p)
mu2 <- mu[p]

# Create the chain
mvn_gibbs <- stableGR:::mvn_gibbs
out.gibbs1 <- mvn_gibbs(N = N, mu = mu, sigma = sigma, p = p)
obj <- list(out.gibbs1)

# Convert to MCMC objects
#out1 <- mcmc(out.gibbs1)
#obj <- mcmc.list(out1)

# Calculate psrfs with gr.diag
results1 <- gr.diag(obj)

##############################################
# Calculate psrfs by for this chain by hand
# Isolate the first component of the  chain
chain1 <- out.gibbs1[ ,1]

# Calculate tau^2 for each chain
tausq1 <- (mcse(chain1, method = "lug")$se)^2 * N 
tausq <- tausq1 
#checked, right

# Calulate s^2 for each chain
sampvar1 <- var(chain1)
ssquared <- sampvar1 
#checked, right

# Calculate sigma^2 estimate
sigsq <- ((N-1) * ssquared + tausq)/N
# checked, right



# Calculate the diagnostic
Rhat <- sigsq / ssquared

that <- sqrt(Rhat)

all.equal(as.numeric(that[1]) , as.numeric(results1$psrf[1]))

