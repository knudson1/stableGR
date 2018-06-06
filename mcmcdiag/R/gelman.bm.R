#####################################
### Function gelman.bm incorporates batch means variance estimates
### into the Gelman-Rubin diagnostic framework 
### x: a list of chains
### the rest of the arguments are leftover from gelman.diag
### Returns the potential scale reduction factor (psrf) 
### for each variable. Note psrf = sqrt(R-hat).
#####################################

gelman.bm <-
function (x, confidence = 0.95, transform = FALSE, df = TRUE,  
    mapping = "determinant", autoburnin = FALSE, 
    multivariate = TRUE, method = "bm", blather = FALSE) 
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

	# Calculate V, the scale of the t distribution.    
	V <- sigsq +  tau2/(Niter * Nchain) 

	df.adj <- 1

    if(Nchain > 1 & df){

		tau2var <- apply(tau2i, 1, var)/Nchain # Calculate the variance of our estimate

		# Calculate cov(tau^2, s^2). Calculate it for each chain, then average.
		# This replaces cov(w,b) from GR. w = s^2 but tau^2 != b.
		cov.s2t2 <- sapply(1:Nvar, function(index) {cov(x = tau2i[index, ], y = s2[index, ])}) 
		# cov(s^2, tau^2) = sample cov(s^2, tau^2) 



		# Variance of the variances, divided by m chains. 
		var.s2 <- apply(s2, 1, var)/Nchain # Known as var.w in GR.


	    # Calculate the variance of V.
		var.V <- ((Niter - 1)/Niter)^2 * var.s2 + 
			+ ((Nchain + 1)/(Nchain * Niter))^2 * tau2var  +  
			+ 2 * (Nchain + 1) * (Niter -1) * cov.s2t2 / (Nchain^2 * Niter^2)

		# Eight, calculate degrees of freedom for our T dist.

		df.V <- (2 * V^2)/var.V 
   		df.adj <- (df.V + 3)/(df.V + 1) 
		# Adjustment keeps SRF finite. 
	}


    #if(df == FALSE) {df.adj <- 1}


	arrr <- V * df.adj / Ssq
	psrf <- sqrt(arrr)

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
        }
        else{
            Sinv <- qr.solve(W)
            thirdpiece <- max(eigen(Sinv %*% Tee, symmetric = FALSE, only.values = TRUE)$values)
        }

		mpsrf <- sqrt(firstpiece + secondpiece*thirdpiece)

	}


	if(blather){
		blather <- list(muhat = muhat, method = method, Niter = Niter, Nchain = Nchain, Nvar = Nvar,
			tausq = tau2, ssq <- s2, sigsq = sigsq, df.adj = df.adj, S = W, Tee = Tee)
	}

	list(psrf = psrf, mpsrf = mpsrf, blather = blather)

   
}


gettau <- function(x1, method, Niter, Nvar) 
{
	asym.var <- numeric(length = Nvar)
	if(method == "bm")
	{
		asym.var <- (mcse.mat(x1, method = method)[ ,2])^2 
	}
	if(method == "wbm")
	{
		bn <- floor(sqrt(Niter))
		asym.var <- 2*(mcse.mat(x1, method = "bm")[ ,2])^2  - (mcse.mat(x1, method = "bm", size = ceiling(bn/2))[ ,2])^2 
	}
	return(asym.var)
}


getT <- function(x, method, Niter, Nvar) 
{

	asym.var <- matrix(0, nrow = Nvar, ncol = Nvar)

	if(method == "bm")
	{
		asym.var <- mcse.multi(x, method = "bm")$cov
	}
	if(method == "wbm")
	{
		bn <- floor(sqrt(Niter))
		asym.var <- 2*mcse.multi(x, method = "bm")$cov - mcse.multi(x, method = "bm", size = ceiling(bn/2))$cov
	}
	return(asym.var)

}

mcse.mat <- mcmcse:::mcse.mat
mcse.mult <- mcmcse:::mcse.multi
gelman.transform <- coda:::gelman.transform

