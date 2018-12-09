#' A Gelman-Rubin diagnostic with batch means
#' 
#' This function uses batch means estimators to calculate a convergence diagnostic for Markov chain Monte Carlo in the spirit of Gelman-Rubin. A univariate `potential scale reduction factor' (PSRF) is calculated for each variable in \code{x}. For multivariate chains, a multivariate PSRF is calculated to take into account the interdependence of the chain's components.  The PSRFs decreases to 1 as the chain length increases. When the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has converged to to the target distribution
#'
#' @param x an \code{mcmc.list} object with more than one chain, and with starting values that are overdispersed with respect to the posterior distribution.
#' @param mapping the function used to map the covariance matrix to a scalar. This is one of \dQuote{\code{determinant}} (determinant of the covariance matrix, the default) or \dQuote{\code{maxeigen}} (the largest eigenvalue of the covariance matrix).
#' @param multivariate a logical flag indicating whether the multivariate potential scale reduction factor should be calculated for multivariate chains
#' @param method the method used to compute the standard error of the chains. This is one of \dQuote{\code{bm}} (batch means, the default), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window)
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#' @param blather a logical flag indicating whether to include additional output
#'
#' @return  \item{psrf}{A vector containing the point estimates of the potential scale reduction factor.}
#' @return \item{mpsrf}{A scalar point estimate of the multivariate potential scale reduction factor.}
#' @return \item{means}{A vector containing the sample means based on the chains provided.}
#' @return \item{blather}{Either \code{FALSE} or a list containing intermediate calculations.}
#'
#' @section Theory: Gelman and Rubin (1992) and Brooks and Gelman (1998) first constructed the univariate and 
#' multivariate potential scale reduction factors (PSRF), respectively, used to diagnose Markov chain 
#' convergence. The function \code{gr.diag} stabilizes the PSRF and improves the PSRF's efficiency by 
#' incorporating batch means estimators for the target variance. The PSRF decreases to 1 as the chain length
#'  increases; when the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has 
#'  converged to to the target distribution. A PSRF convergence threshold can be calculated using 
#'  \code{choosepsrf}.
#'
#'
#' @section References:
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. \cr
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. \cr
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @export
#'

gr.diag <-
function (x, mapping = "determinant",  
    multivariate = TRUE, method = "lug", autoburnin = FALSE, blather = FALSE) 
{
    x <- as.mcmc.list(x)
    
    if (autoburnin && start(x) < end(x)/2) 
      x <- window(x, start = end(x)/2 + 1)

	# Define some notation.
    Niter <- niter(x)  # number of iterations per chains. We also call this n.
    Nchain <- nchain(x) # number of chains. We also call this m.
    Nvar <- nvar(x) # number of variables
    xnames <- varnames(x)

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
    W <- apply(Si2, c(1, 2), mean) # Average the vcov matrices across chains. #aka S
    Ssq <- diag(W) # Isolate the variances, throw away covariances.
	# For each chain, find sample variance for each variable. 
	s2 <- matrix(apply(Si2, 3, diag), nrow = Nvar, ncol = Nchain)

	# Third, calculate tau^2 and its variance for each variable. 
	# This replaces the GR b.
	# Sample variance of the sample means (between chain vars) calculated using batch means.
  tau2i <- matrix(sapply(x, gettau, method = method)*Niter,  ncol = Nchain) # For each chain
	tau2 <- apply(tau2i, 1, mean)  # Average over the chains

	# Calculate the estimate of sigma^2.
	sigsq <- (Niter - 1) * Ssq/Niter + tau2 / Niter 

	arrr <- sigsq / Ssq
	psrf <- sqrt(arrr)
	
	blatherout <- blather
	
	if(blather){
	  blatherout <- list(muhat = muhat, method = method, Niter = Niter, Nchain = Nchain, Nvar = Nvar,
	                  tausq = tau2, ssq <- s2, sigsq = sigsq) }

	mpsrf <- multivariate
	
	if(multivariate && Nvar > 1){
		Ti <- lapply(x, getT, method = method)  # For each chain
		Tee <- matrix(Reduce("+", Ti)  / Nchain, nrow = Nvar)
		Tee <- adjust_matrix(Tee, Niter)

		firstpiece <- (Niter-1)/Niter
		secondpiece <- 1/Niter
		mango <- solve(W, Tee) #S^{-1}T
		eigs <- eigen(mango, symmetric = FALSE, only.values = TRUE)$values
		
		thirdpiecedet <- (prod(eigs))^(1/Nvar)
		thirdpiecemax <- max(eigs)
		
		mpsrfdet <- sqrt(firstpiece + secondpiece*thirdpiecedet)
		mpsrfmax <- sqrt(firstpiece + secondpiece*thirdpiecemax)		
		

        if(mapping == "determinant"){ mpsrf <- mpsrfdet
        }else{ mpsrf <- mpsrfmax
        }


		if(blather){
		  blatherout$S <- W 
		  blatherout$Tee <- Tee
		  blatherout$eigenvalues <- eigs
		  blatherout$mpsrfdet <- mpsrfdet
		  blatherout$mpsrfmax <- mpsrfmax
		}
	}




	list(psrf = psrf, mpsrf = mpsrf, means = muhat, blather = blatherout)

   
}


gettau <- function(x1, method) 
{
	(mcse.mat(x1, method = method)[ ,2])^2 

}


getT <- function(x, method) 
{
  mcse.multi(x, method = method)$cov
  # foo <- mcse.multi(x, method = method)
	# Tee <- adjust_matrix(foo$cov, N = foo$nsim)
	# Tee
}

adjust_matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 1/2)
{
  mat.adj <- mat
  adj <- epsilon*N^(-b)
  #for(i in 1:ncol(mat.adj)){mat.adj[i,i] <- }
  diag(mat.adj) <- pmax(adj, diag(mat.adj))
  vars <- diag(mat.adj)
  corr <- cov2cor(mat.adj)
  eig <- eigen(corr)
  adj.eigs <- pmax(eig$values, adj)
  mat.adj <- diag(vars^(1/2))%*% eig$vectors %*% diag(adj.eigs) %*% t(eig$vectors) %*% diag(vars^(1/2))
  return(mat.adj)
}

mcse.mat <- mcmcse:::mcse.mat
mcse.multi <- mcmcse:::mcse.multi
gelman.transform <- coda:::gelman.transform

