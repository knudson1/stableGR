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

withfunc <- stable.GR(chainlist, method = "lug", blather = TRUE)
blatherout <- withfunc$blather

# calculate batch size
bees <- sapply(chainlist, batchSize, method = "bm")
b <- floor(mean(bees))
a <- floor(N/b)
Nneeded <- a*b
Ngoaway <- N - Nneeded

all.equal(b, blatherout$b)
all.equal(a, blatherout$a)
all.equal(Nneeded, blatherout$Nneeded)

#trim away beginning of each chain
if(Ngoaway ==0) trimmedchains <- chainlist
if(Ngoaway >0) {
    trimmedchain1 <- out.gibbs1[-(1:Ngoaway),]
    trimmedchain2 <- out.gibbs2[-(1:Ngoaway),]
    trimmedchains <- list(trimmedchain1, trimmedchain2)
}

#stack the trimmed chains
stacked <- do.call(rbind, trimmedchains)
all.equal(stacked, blatherout$stackedchains)

# Calculate tau^2 by hand but using mcse.mat
mat <- mcse.mat(stacked, method = "lug", size = b)
#SE are in second column, and we want first comp
tausq <- (mat[,2])^2 * Nneeded
all.equal(tausq, blatherout$tausq)


# Calulate s^2 (variance for each component) without apply
ssquared <- rep(0, ncol(stacked))
for(i in 1:ncol(stacked)){
    ssquared[i] <- var(stacked[,i])
}


# Calculate sigma^2 estimate
sigsq <- ((Nneeded-1) * ssquared + tausq)/(Nneeded)
all.equal(sigsq, blatherout$sigmasq)

# Calculate the diagnostic
Rhat <- sigsq / ssquared
psrfbyhand <- sqrt(Rhat)
names(psrfbyhand) <- NULL

#moment of truth
all.equal(as.numeric(withfunc$psrf), as.numeric(psrfbyhand))



