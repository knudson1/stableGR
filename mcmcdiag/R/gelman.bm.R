#####################################
### Incorporates lugsail batch means variance 
### estimatesinto the Gelman-Rubin diagnostic framework 
### x: a list of chains
### the rest of the arguments are leftover from gelman.diag
### Returns the potential scale reduction factor (psrf) 
### for each variable. Note psrf = sqrt(R-hat).
#####################################

gelman.bm <-
function (x, confidence = 0.95, transform = FALSE,  
    mapping = "determinant", autoburnin = FALSE, 
    multivariate = TRUE, method = "lug", blather = FALSE) 
{
    x <- as.mcmc.list(x)
    if (autoburnin && start(x) < end(x)/2) 
        x <- window(x, start = end(x)/2 + 1)

	# Define some notation.
    Niter <- niter(x)  # number of iterations per chains. We also call this n.
    Nchain <- nchain(x) # number of chains. We also call this m.
    Nvar <- nvar(x) # number of variables
    xnames <- varnames(x)

	# Transform to logit or log if asked and if applicable.
    if (transform) 
        x <- gelman.transform(x) 

	# Since x is a list of markov chains, turn each into matrix.
    x <- lapply(x, as.matrix) 

	# First, calculate sample means.
	# Calculate column means for each chain.
    xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE),nrow = Nvar, ncol = Nchain) 
	# Average the means across chains
    muhat <- apply(xbar, 1, mean) 

	# Second, calculate overall sample variance for each variable.
	# Calculate vcov matrix for variables in each chain. List length = nchain.
    Si2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, Nvar, Nchain)) 
    W <- apply(Si2, c(1, 2), mean) # Average the vcov matrices across chains.
    Ssq <- diag(W) # Isolate the variances, throw away covariances.
	# For each chain, find sample variance for each variable. 
	s2 <- matrix(apply(Si2, 3, diag), nrow = Nvar, ncol = Nchain)

	# Third, calculate tau^2 and its variance for each variable. 
	# This replaces the GR b.
	# Sample variance of the sample means (between chain vars) calculated using batch means.
    tau2i <- matrix(sapply(x, gettau, method = method, Niter = Niter, Nvar = Nvar)*Niter,  ncol = Nchain) # For each chain
	tau2 <- apply(tau2i, 1, mean)  # Average over the chains

	# Calculate the estimate of sigma^2.
	sigsq <- (Niter - 1) * Ssq/Niter + tau2 / Niter 

	arrr <- sigsq / Ssq
	psrf <- sqrt(arrr)
	
	if(blather){
	  blatherout <- list(muhat = muhat, method = method, Niter = Niter, Nchain = Nchain, Nvar = Nvar,
	                  tausq = tau2, ssq <- s2, sigsq = sigsq) }

	if(multivariate == TRUE  && Nvar > 1){
		Ti <- lapply(x, getT, method = method, Niter = Niter, Nvar = Nvar)  # For each chain
		Tee <- matrix(Reduce("+", Ti)  / Nchain, nrow = Nvar)

		firstpiece <- (Niter-1)/Niter
		secondpiece <- (Nchain+1)/(Nchain*Niter)

        if(mapping == "determinant"){
    		eigenS <- eigen(W, symmetric=TRUE, only.values = TRUE)$values
    		bottom <- exp(mean(log(eigenS)))
    
    		eigenT <- eigen(Tee, symmetric = TRUE, only.values = TRUE)$values
    		top <- (prod(eigenT))^(1/Nvar)
    
    		thirdpiece <- (top/bottom)
        }else{
            Sinv <- qr.solve(W)
            thirdpiece <- max(eigen(Sinv %*% Tee, symmetric = FALSE, only.values = TRUE)$values)
        }

		mpsrf <- sqrt(firstpiece + secondpiece*thirdpiece)

	}


	if(blather){
		blatherout$S <- W 
		blatherout$Tee <- Tee
		blather <- blatherout
	}

	list(psrf = psrf, mpsrf = mpsrf, means = muhat, blather = blather)

   
}


gettau <- function(x1, method, Niter=NULL, Nvar) 
{
	asym.var <- numeric(length = Nvar)
	asym.var <- (mcse.mat(x1, method = method)[ ,2])^2 
# 	if(method == "bm")
# 	{
# 		asym.var <- (mcse.mat(x1, method = method)[ ,2])^2 
# 	}
# 	if(method == "wbm")
# 	{
# 		bn <- floor(sqrt(Niter))
# 		asym.var <- 2*(mcse.mat(x1, method = "bm")[ ,2])^2  - (mcse.mat(x1, method = "bm", size = ceiling(bn/2))[ ,2])^2 
# 	}
	return(asym.var)
}


getT <- function(x, method, Niter, Nvar) 
{
	#asym.var <- matrix(0, nrow = Nvar, ncol = Nvar)
	foo <- mcse.multi(x, method = method)
	Tee <- adjust_matrix(foo$cov, N = foo$nsim)
	Tee
}

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 1/2)
{
  mat.adj <- mat
  adj <- epsilon*N^(-b)
  vars <- diag(mat)
  corr <- cov2cor(mat)
  eig <- eigen(corr)
  adj.eigs <- pmax(eig$values, adj)
  mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
  return(mat.adj)
}

mcse.mat <- mcmcse:::mcse.mat
mcse.mult <- mcmcse:::mcse.multi
gelman.transform <- coda:::gelman.transform

