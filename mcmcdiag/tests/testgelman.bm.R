set.seed(10)
library(mvtnorm)
library(mcmcse)
library(coda)
library(mcmcdiag)

################ 
# Start by making a few chains to work with

# Details on the chain construction
p <- 5
N <- 100
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
# Just for the first variable

# Write the psrf for the specific obj chains
# Isolate the first component of the two chains
chain1 <- out1[ ,1]
chain2 <- out2[ ,1]

# Calculate tau^2 for each chain
tausq1 <- (mcse(chain1)$se)^2 * N 
tausq2 <- (mcse(chain2)$se)^2 * N
(tausq <- .5 * (tausq1 + tausq2))
#checked, right

# Calulate s^2 for each chain
sampvar1 <- var(chain1)
sampvar2 <- var(chain2)
(ssquared <- .5 * (sampvar1 + sampvar2))
#checked, right

# Calculate sigma^2 estimate
sigsq <- ((N-1) * ssquared + tausq)/N
#checked, right

# Calculate Vhat 
Vhat <- sigsq + tausq/(N * 2)
#checked, right

# Moving from point estimates to the vars of estimators
vartausq <- var(c(tausq1, tausq2))/2
varssq <- var(c(sampvar1, sampvar2))/2
thecov <- cov(x=c(sampvar1, sampvar2), y=c(tausq1, tausq2))
thingy1 <- (sampvar1 - ssquared)*(tausq1 - tausq)
thingy2 <- (sampvar2 - ssquared)*(tausq2 - tausq)
altcov <- (sampvar1 - ssquared)*(tausq1 - tausq) + (sampvar2 - ssquared)*(tausq2 - tausq)
all.equal(altcov, thecov) #just to be sure!!

# Variance of Vhat
m <- 2
varv <- ((N-1)/N)^2 * varssq + ((m+1)/(m*N))^2 * vartausq + 2*(m+1)*(N-1)*thecov/(m^2 * N^2)
# checked, looks great

# Calculate df for the T distribution
defo <- 2 * Vhat^2 / varv
dfchunk <- (defo + 3)/(defo + 1)
# checked


# Calculate the diagnostic
Rhat <- Vhat * dfchunk / ssquared

that <- sqrt(Rhat)

byhand <- list(tausq1 = tausq1, tausq2 = tausq2, tausq = tausq, sampvar1 = sampvar1, ssquared = ssquared, sigsq = sigsq, Vhat = Vhat, vartausq = vartausq, varssq = varssq, thecov = thecov, varv = varv, defo = defo, Rhat = Rhat)

byhand



(withfunc <- gelman.bm(obj))
withfunc$psrf[1]
all.equal(as.numeric(withfunc$psrf[1]), that)



