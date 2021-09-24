#' Asymptotic covariance matrix estimation for Markov chain Monte Carlo
#' 
#' Estimates the asymptotic covariance matrix for Monte Carlo estimators,  compatible with multiple chains. If a single chain is input, it calls \code{mcmcse::mcse.multi}. 
#' 
#' @param x a list of matrices, where each matrix is \eqn{n \times p}. Each row of the matrices represents one step of the chain. Each column of the matrices represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object.
#' @param multivariate a logical flag indicating whether the full matrix is returned (TRUE) or only the diagonals (FALSE)
#' @param method the method used to compute the matrix. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size options are \code{NULL} (default, which calculates an ideal batch size), character values of \code{sqroot} and \code{cuberoot}, or any numeric value between 1 and \eqn{n}. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#' @param adjust this argument is now obselete due to package updates.
#'
#' @return  The asymptotic variance estimate (if \code{multivariate = FALSE}) or the asymptotic covariance matrix (if \code{multivariate = TRUE}) in the Markov chain central limit theorem. 
#'
#' @details The function returns estimate of the univariate or multivariate asymptotic (co)variance of Monte Carlo estimators. If \eqn{X_1, \dots X_n} are the MCMC samples, then function returns the estimate of \eqn{\lim_{n\to \infty} n Var(\bar{X})}. In other words, if a Markov chain central limit holds such that, as \eqn{n \to \infty}
#'   \deqn{\sqrt{n}(\bar{X} - \mu) \to N(0, \Sigma) }
#' then the function returns an estimator of \eqn{\Sigma} from the m different chains. If \code{multivariate == FALSE}, then only the diagonal of \eqn{\Sigma} are returned.
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
#' asym.var(x)
#' 
#' @export
#'


asym.var <- function (x, 
                      multivariate = TRUE, 
                      method = "lug", 
                      size = NULL, 
                      autoburnin = FALSE, 
                      adjust = TRUE) 
{
  # perform various checks on markov chains
  x <- mcmcchecks(x, autoburnin = autoburnin)
  
  #trim the chains to prep for replicated bm
  out <- size.and.trim(x = x, size = size)
  b <- out$b
  Nneeded <- out$Nneeded
  trimmedchains <- out$trimmedchains
  Nvar <- out$Nvar
    
  ## stack the chains into a single matrix
  stackedchains <- do.call(rbind, trimmedchains)
  
  ## calculate tau squared using replicated batch means
  if(multivariate == FALSE){
      Tee <- gettau(stackedchains, method = method, size = b) * Nneeded 
      #### Doots, mcse breaks if b = 1 (says is too small)
  }
  
  ## calculate T using replicated batch means
  if(multivariate){
      if(Nvar == 1)stop("The option multivariate = TRUE requires a Markov chain with multiple variables. If you have a univariate Markov chain, use multivariate = FALSE.")
      Tee <- getT(stackedchains, method = method, size = b)
  }

  Tee

}




