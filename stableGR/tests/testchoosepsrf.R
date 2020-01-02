library(stableGR)
minESS <- mcmcse::minESS
p <- 1
alpha <- .05
m <- 5


#### Check the case where delta is null (epsilon is instead given)
# That is, find the effective sample size based on epsilon
epsilon <- .05
out1a <- target.psrf(p, epsilon = epsilon, m = m, alpha = alpha)
Tee <- as.numeric(minESS(p, epsilon, alpha = alpha)) #min effective sample size
del <- sqrt(1 + m/Tee) - 1
arr <- 1 + del  
out1b <- list(psrf = arr, epsilon = epsilon)
all.equal(out1a, out1b)

#### Check the case where delta is NOT null (it's given)
# That is, find the effective sample size based on delta
delta <- .01
out2a <- target.psrf(p, delta = delta, m = m, alpha = alpha)
Tee <- M <- m/((1+delta)^2 - 1)
epsilon <- as.numeric(minESS(p, eps = epsilon, ess = M))
del <- sqrt(1 + m/Tee) - 1
arr <- 1 + del  
out2b <- list(psrf = arr, epsilon = epsilon)
all.equal(out2a, out2b)

