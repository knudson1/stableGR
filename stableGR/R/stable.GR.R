#' Gelman-Rubin diagnostic using stable variance estimators
#' 
#' This function uses fast and strongly consistent estimators estimators of Monte Carlo variance to calculate the Gelman-Rubin convergence diagnostic for Markov chain Monte Carlo. A univariate `potential scale reduction factor' (PSRF) is calculated for each variable in \code{x}. For multivariate chains, a multivariate PSRF is calculated to take into account the interdependence of the chain's components.  The PSRFs decrease to 1 as the chain length increases. When the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has converged to the target distribution.
#'
#' @param x a list of matrices, where each matrix represents one Markov chain sample. Each row of the matrices represents one step of the chain. Each column of the matrices represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object.
#' @param multivariate a logical flag indicating whether the multivariate potential scale reduction factor should be calculated for multivariate chains.
#' @param mapping the function used to map the covariance matrix to a scalar. This is one of \dQuote{\code{determinant}} (determinant of the covariance matrix, the default) or \dQuote{\code{maxeigen}} (the largest eigenvalue of the covariance matrix).
#' @param method the method used to compute the standard error of the chains. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size options are \code{NULL} (default, which calculates an ideal batch size), character values of \code{sqroot} and \code{cuberoot}, or any numeric value between 1 and \eqn{n}. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#' @param blather a logical flag indicating whether to include additional output.
#'
#' @return \item{psrf}{A vector containing the point estimates of the PSRF.}
#' @return \item{mpsrf}{A scalar point estimate of the multivariate PSRF.}
#' @return \item{means}{A vector containing the sample means based on the chains provided.}
#' @return \item{n.eff}{A scalar point estimate of the effective sample size.}
#' @return \item{blather}{Either \code{FALSE} or a list containing intermediate calculations.}
#'
#' @section Theory: Gelman and Rubin (1992) and Brooks and Gelman (1998) first constructed the univariate and 
#' multivariate potential scale reduction factors (PSRF), respectively,  to diagnose Markov chain 
#' convergence. The function \code{stable.GR} stabilizes the PSRF and improves the PSRF's efficiency by 
#' incorporating lugsail estimators for the target variance. The PSRF decreases to 1 as the chain length
#'  increases; when the PSRF becomes sufficiently close to 1, the sample collected by the Markov chain has 
#'  converged to to the target distribution. A PSRF convergence threshold can be calculated using 
#'  \code{choosepsrf}.
#'
#' @examples 
#' library(stableGR)
#' set.seed(100)
#' p <- 2
#' n <- 100 # n is tiny here purely for demo purposes.
#' # use n much larger for real problems!
#' 
#' sig.mat = matrix(c(1, .8, .8, 1), ncol = 2, nrow = 2)
#'
#' # Making 3 chains
#' chain1 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain2 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain3 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#'
#' # find GR diagnostic using all three chains
#' x <- list(chain1, chain2, chain3)
#' stable.GR(x) 
#' 
#'
#' @section References:
#' Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384. 
#' 
#' Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.
#' 
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. \cr
#'
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. \cr
#'
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @export
#'

stable.GR <-
function (x, 
          multivariate = TRUE, 
          mapping = "determinant",  
          method = "lug", 
          size = NULL, 
          autoburnin = FALSE, 
          blather = FALSE) 
{
  # make sure markov chains pass various checks
  x <- mcmcchecks(x, autoburnin = autoburnin)
  xnames <- colnames(x[[1]])
  Nchain <- length(x)
  
  # preliminary number crunching
  out <- size.and.trim(x=x, size = size)
  Nvar <- out$Nvar
  Nneeded <- out$Nneeded
  a <- out$a
  b <- out$b
  trimmedchains <- out$trimmedchains
  # chains prepared for rep bm, still in list
  
  #stack the trimmed chains so it looks like one chain and the batches will line up for rep bm
  stackedchains <- do.call(rbind, trimmedchains)
  
  
	# First, calculate sample means for each component.
  muhat <- apply(stackedchains, 2, mean)

	# Second, calculate overall sample variance for each variable.
	# Calculate vcov matrix for variables in each chain. List length = nchain.
  W <- var(stackedchains)
  Ssq <- diag(W) # Isolate the variances, throw away covariances.

	# Third, calculate tau^2, the sample variance of the sample means (between chain vars)
  tau2 <- asym.var(x, multivariate = FALSE, method = method, size = size, autoburnin = FALSE)
	
	# Calculate the estimate of sigma^2.
	sigsq <- (Nneeded - 1) * Ssq/Nneeded + tau2 / Nneeded 

	arrr <- sigsq / Ssq
	psrf <- sqrt(arrr)

	blatherout <- blather
	
	if(blather){
	  blatherout <- list(method = method, 
	                     Nneeded = Nneeded,
	                     Nchain = Nchain, 
	                     Nvar = Nvar,
	                     asymVars = tau2, 
	                     sigmasq = sigsq,
	                     a = a,
	                     b = b,
	                     stackedchains = stackedchains,
	                     tausq = tau2) }

	denom <- arrr - ((Nneeded-1)/Nneeded)
	n.eff <- Nchain/denom
	
	mpsrf <- multivariate
	
	if(multivariate && Nvar > 1){
	  Tee <- asym.var(x, multivariate = TRUE, method = method, size = size, autoburnin = FALSE)

		firstpiece <- (Nneeded-1)/Nneeded
		secondpiece <- 1/Nneeded
		mango <- solve(W, Tee) #S^{-1}T
		eigs <- eigen(mango, symmetric = FALSE, only.values = TRUE)$values
		
		thirdpiecedet <- (prod(eigs))^(1/Nvar)
		thirdpiecemax <- max(eigs)
		
		mpsrfdet <- sqrt(firstpiece + secondpiece*thirdpiecedet)
		mpsrfmax <- sqrt(firstpiece + secondpiece*thirdpiecemax)		
		

        if(mapping == "determinant"){ mpsrf <- mpsrfdet
        }else{ mpsrf <- mpsrfmax
        }

		denom <- mpsrf^2 - ((Nneeded-1)/Nneeded)
		n.eff <- Nchain/denom

		if(blather){
		  blatherout$S <- W 
		  blatherout$AsymVarMatrix <- Tee
		  blatherout$eigenvalues <- eigs
		  blatherout$mpsrfDet <- mpsrfdet
		  blatherout$mpsrfMaxEigen <- mpsrfmax
		}
	}

  names(muhat) <- xnames

	list(psrf = psrf, mpsrf = mpsrf, means = muhat, n.eff = n.eff, blather = blatherout)

   
}


