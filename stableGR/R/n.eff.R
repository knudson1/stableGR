#' Effective sample size 
#' 
#' For an estimator, effective sample size is the number of independent samples with the same standard error as a correlated sample. This function calculates effective sample size for a set of Markov chains using lugsail variance estimators. This also determines whether the Markov chains have converged. If they have not, this function approximates the number of samples needed.
#'
#' @param x a list of matrices, where each matrix represents one Markov chain sample. Each row of the matrices represents one step of the chain. Each column of the matrices represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object.
#' @param multivariate a logical flag indicating whether the effective sample size should be calculated for multivariate chains.
#' @param epsilon relative precision level. Values less than .10 are recommended.
#' @param delta desired delta value - the cutoff for potential scale reduction factor. 
#' @param alpha significance level for confidence regions for the Monte Carlo estimators.
#' @param method the method used to compute the standard error of the chains. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size can take character values of \dQuote{sqroot} and \dQuote{cuberoot} or any numeric value between 1 and n. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.

#'
#' @return \item{n.eff}{a scalar point estimate of the effective sample size.}
#' @return \item{converged}{a logical indicating whether sufficient samples have been obtained.}
#' @return \item{n.target}{NULL (if \code{converged == TRUE}) or a scalar estimate of the chain length required for convergence, assuming the number of chains is unchanged.  }
#'
#' @examples 
#' library(stableGR)
#' set.seed(100)
#' p <- 2
#' n <- 10000
#'
#' sig.mat = matrix(c(1, .8, .8, 1), ncol = 2, nrow = 2)
#' # Making 3 chains
#' \dontrun{ chain1 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain2 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain3 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#'
#' # find ESS using all three chains
#' x <- list(chain1, chain2, chain3)
#' n.eff(x) 
#' }

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


n.eff <- function(x,  multivariate = TRUE, epsilon = .05, delta = NULL, alpha = .05, method = "lug", size = "sqroot", 
                  autoburnin = FALSE){ 
  
  # make sure markov chains pass various checks
  x <- mcmcchecks(x, autoburnin = autoburnin)
  
  # Define some notation.
  Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
  Nvar <- ncol(x[[1]]) # number of variables
  Nchain <- length(x)
  
  #univariate case
  W <- s.hat(x)
  Ssq <- diag(W) 
  tau2 <- asym.var(x, method = method, size = size, autoburnin = FALSE)
  currentESS <- Nchain * Niter * (Ssq/tau2)^(1/Nvar)  
  
    
  if(multivariate && Nvar > 1){
    Tee <- asym.var.mat(x, method = method, size = size, autoburnin = FALSE, adjust = TRUE)
    mango <- solve(Tee, W) #S T^{-1}
    eigs <- eigen(mango, symmetric = FALSE, only.values = TRUE)$values
    detpiece <- (prod(eigs))^(1/Nvar)
    currentESS <- Nchain * Niter * detpiece
    }
  
  # prepare and do comparison to our goal
  p <- Nvar
  m <- Nchain
  targ <- target.psrf(p = p, m = m, epsilon = epsilon, delta = delta, alpha = alpha)
  targetESS <- RtoESS(targ$psrf, m)
  converged <- FALSE
  if(currentESS >= targetESS){ converged <- TRUE}
  
  # if the sample hasn't converged, approximate the sample size that would result in convergence
  ntarget <- nmore <- NULL
  if(converged == FALSE){
    ncurrent <- nrow(x[[1]])
    ntarget <- ceiling(ncurrent*targetESS/currentESS)
    nmore <- ceiling(ntarget - ncurrent)
  }  
  
  list(n.eff = currentESS, converged = converged, n.target = ntarget)
}

RtoESS <- function(arrr,  m){
  m/(arrr^2 - 1)
}


