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

################ 
# Perform unit test using the two chains 
# Just for the first variable

chainlist <- list(out.gibbs1, out.gibbs2)

# calculate batch size
bees <- sapply(chainlist, batchSize, method = "bm")
b <- floor(mean(bees))
a <- floor(N/b)
Nneeded <- a*b
Ngoaway <- N - Nneeded

#trim away beginning of each chain
if(Ngoaway ==0) trimmedchains <- chainlist
if(Ngoaway >0) {
    trimmedchain1 <- out.gibbs1[-(1:Ngoaway),]
    trimmedchain2 <- out.gibbs2[-(1:Ngoaway),]
    trimmedchains <- list(trimmedchain1, trimmedchain2)
}

#stack the trimmed chains
stacked <- do.call(rbind, trimmedchains)

# Calculate tau^2 for first component by hand
mat <- mcse.mat(stacked, method = "lug", size = b)
#SE are in second column, and we want first comp
tausq <- (mat[1,2])^2 * 2 * Nneeded
names(tausq) <- c('se')

#calculate tau^2 for first component using asym.var
onecomp <- stacked[,1]
onecomp <-matrix(onecomp, ncol=1)
coffee <- asym.var(onecomp, multivariate = FALSE, method = "lug", size = b)

#do two tau^2 calculations match?
all.equal(coffee, tausq)

# Calulate s^2 (combine all samples, calc var)
ssquared <- var(onecomp)

# Calculate sigma^2 estimate
sigsq <- ((Nneeded-1) * ssquared + tausq)/(Nneeded)



# Calculate the diagnostic
Rhat <- sigsq / ssquared

that <- sqrt(Rhat)
names(that) <- NULL



withfunc <- stable.GR(chainlist, method = "lug")
all.equal(as.numeric(withfunc$psrf[1]), as.numeric(that))



