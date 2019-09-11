#' Asymptotic covariance matrix estimation for Markov chain Monte Carlo
#' 
#' Estimates the asymptotic covariance matrix for Monte Carlo estimators,  compatible with multiple chains. If a single chain is input, it calls \code{mcmcse::mcse.multi}. 
#' 
#' @param x a list of matrices, where each matrix is \eqn{n \times p}. Each row of the matrices represents one step of the chain. Each column of the matrices represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object.
#' @param multivariate a logical flag indicating whether the full matrix is returned (TRUE) or only the diagonals (FALSE)
#' @param method the method used to compute the matrix. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size can take character values of \code{sqroot} and \code{cuberoot} or any numeric value between 1 and \eqn{n}. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#' @param adjust a logical flag indicating whether the covariance matrix should be adjusted, when necessary, to ensure it is positive-definite.
#'
#' @return  The asymptotic variance estimate (if \code{multivariate = FALSE}) or the asymptotic covariance matrix (if \code{multivariate = TRUE}) in the Markov chain central limit theorem. 
#'
#' @details The function returns estimate of the univariate or multivariate asymptotic (co)variance of Monte Carlo estimators. If \eqn{X_1, \dots X_n} are the MCMC samples, then function returns the estimate of \eqn{\lim_{n\to \infty} n Var(\bar{X})}. In other words, if a Markov chain central limit holds such that, as \eqn{n \to \infty}
#'   \deqn{\sqrt{n}(\bar{X} - \mu) \to N(0, \Sigma) }
#' then the function returns an estimator of \eqn{\Sigma} from the m different chains. If \code{multivariate == FALSE}, then only the diaongal of \eqn{\Sigma} are returned.
#' @section References:
#' Vats, D. and Knudson, C. Revisiting the Gelman-Rubin Diagnostic.	arXiv:1812.09384. 
#' 
#' Vats, D. and Flegal, J. Lugsail lag windows and their application to MCMC. arXiv: 1809.04541.
#' 
#' Flegal, J. M. and Jones, G. L. (2010) Batch means and spectral variance estimators in Markov chain Monte Carlo. \emph{The Annals of Statistics}, \bold{38}, 1034--1070.
#'
#' Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, \emph{Statistical Science}, \bold{7}, 457-511. 
#'
#' Brooks, SP. and Gelman, A. (1998) General methods for monitoring convergence of iterative simulations. \emph{Journal of Computational and Graphical Statistics}, \bold{7}, 434-455.
#'
#' @examples 
#' library(stableGR)
#' set.seed(100)
#' p <- 2
#' n <- 1000
#'
#' sig.mat = matrix(c(1, .8, .8, 1), ncol = 2, nrow = 2)
#'
#' # Making 3 chains
#' \dontrun{chain1 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain2 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#' chain3 <- mvn.gibbs(N = n, p = p, mu = rep(1,p), sigma = sig.mat)
#'
#' # find GR diagnostic using all three chains
#' x <- list(chain1, chain2, chain3)
#' asym.var(x)
#' }
#' @export
#'


asym.var <- function (x, multivariate = TRUE, method = "lug", size = "sqroot", autoburnin = FALSE, adjust = TRUE) 
{
  mcse.mat <- mcmcse::mcse.mat
  mcse.multi <- mcmcse::mcse.multi
  
  # perform various checks on markov chains
  x <- mcmcchecks(x, autoburnin = autoburnin)
  
  # Define some notation
  Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
  Nvar <- ncol(x[[1]]) # number of variables
  Nchain <- length(x)
  
  # When we have multiple chains, we need to do replicated batch means
  # meaning we need to calculate the batch sizes manually
  if (size == "sqroot") {
      b = floor(sqrt(Niter))
      a = floor(Niter/b)
  }  else if (size == "cuberoot") {
      b = floor(Niter^(1/3))
      a = floor(Niter/b)
  }  else {
      if (!is.numeric(size) || size <= 1 || size == Inf) 
          stop("'size' must be a finite numeric quantity larger than 1.")
      b = floor(size)
      a = floor(Niter/b)
  }

  
  ## trim away beginnings of each chain (if necessary)
  Nneeded <- a*b
  Ntrim <- Niter - Nneeded
  trimmedchain <- x
  if(Ntrim > 0){
      for(i in 1:Nchain){
          removethese <- 1:Ntrim
          trimmedchain[[i]] <- matrix(x[[i]][-removethese,], ncol = Nvar)
      }
  }
  
  ## stack the chains into a single matrix
  stackedchains <- do.call(rbind,trimmedchain)
  
  ## calculate tau squared using replicated batch means
  if(multivariate == FALSE){
      #tau2i <- matrix(sapply(stackedchains, gettau, method = method, size = b)*Niter, ncol = Nchain)
      #tau2 <- apply(tau2i, 1, mean)
      Tee <- gettau(stackedchains, method = method, size = b) * Niter * Nchain
  }
  
  ## calculate T using replicated batch means
  if(multivariate){
      if(Nvar == 1)stop("The option multivariate = TRUE requires a Markov chain with multiple variables. If you have a univariate Markov chain, use multivariate = FALSE.")
      Tee <- getT(stackedchains, method = method, size = b)
      if(adjust == TRUE) Tee <- adjust.matrix(Tee, Niter)      
  }
  
 
  # ## calculate T using old method
  # if(multivariate){
  #     if(Nvar == 1)stop("The option multivariate = TRUE requires a Markov chain with multiple variables. If you have a univariate Markov chain, use multivariate = FALSE.")
  #     
  #     Ti <- lapply(x, getT, method = method, size = b)  # For each chain
  #     Tee <- matrix(Reduce("+", Ti)  / Nchain, nrow = Nvar)
  #     if(adjust == TRUE) Tee <- adjust.matrix(Tee, Niter)      
  # }

  
  Tee

}

