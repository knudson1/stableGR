#' Asymptotic variance estimation for Markov chain Monte Carlo
#' 
#' This function calculates asymptotic variance estimates for Markov chain Monte Carlo.  
#' 
#' @param x a list of matrices, where each matrix represents one Markov chain. Each row represents one step of the chain. Each column represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object. The starting values of the chain(s) should be overdispersed with respect to the posterior distribution.
#' @param method the method used to compute the standard error of the chains. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size can take character values of \dQuote{sqroot} and \dQuote{cuberoot} or any numeric value between 1 and n. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#'
#' @return  The asymptotic variance estimates.
#'
#'
#' @section References:
#' Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384. 
#' 
#' Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.
#' 
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070. \cr
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. \cr
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @export
#'

asym.var <- function (x, method = "lug", size = "sqroot", autoburnin = FALSE) 
{
  # perform various checks on markov chains
  x <- mcmcchecks(x, autoburnin = autoburnin)
  
  # Define some notation
  Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
  Nvar <- ncol(x[[1]]) # number of variables
  Nchain <- length(x)

	# Sample variance of the sample means (between chain vars) 
  tau2i <- matrix(sapply(x, gettau, method = method)*Niter,  ncol = Nchain) # For each chain
	tau2 <- apply(tau2i, 1, mean)  # Average over the chains

	tau2
   
}

mcse.mat <- mcmcse:::mcse.mat
mcse.multi <- mcmcse:::mcse.multi

mcmcchecks <- stableGR:::mcmcchecks
gettau <- stableGR:::gettau
