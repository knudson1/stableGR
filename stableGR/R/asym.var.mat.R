#' Asymptotic variance matrix estimation for Markov chain Monte Carlo
#' 
#' This function estimates the asymptotic variance matrix for Markov chain Monte Carlo.  
#' 
#' @param x a list of matrices, where each matrix represents one Markov chain. Each row represents one step of the chain. Each column represents one variable. A list with a single matrix (chain) is allowed. Optionally, this can be an \code{mcmclist} object. The starting values of the chain(s) should be overdispersed with respect to the posterior distribution.
#' @param method the method used to compute the standard error of the chains. This is one of \dQuote{\code{lug}} (lugsail, the default), \dQuote{\code{bm}} (batch means), \dQuote{\code{obm}} (overlapping batch means), \dQuote{\code{tukey}} (spectral variance method with a Tukey-Hanning window), or \dQuote{\code{bartlett}} (spectral variance method with a Bartlett window).
#' @param size can take character values of \dQuote{sqroot} and \dQuote{cuberoot} or any numeric value between 1 and n. Size represents the batch size in \dQuote{\code{bm}} (batch means) and the truncation point in \dQuote{\code{bartlett}} and \dQuote{\code{tukey}}. sqroot means size is floor(n^(1/2) and cuberoot means size is floor(n^(1/3)).
#' @param autoburnin a logical flag indicating whether only the second half of the series should be used in the computation.  If set to TRUE and \code{start(x)} is less than \code{end(x)/2} then start of series will be adjusted so that only second half of series is used.
#' @param adjust a logical flag indicating whether the covariance matrix should be adjusted, when necessary, to ensure it is positive-definite.
#'
#' @return  The asymptotic variance matrix.
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

asym.var.mat <- function (x, method = "lug", size = "sqroot", autoburnin = FALSE, adjust = TRUE) 
{
  # perform various checks on markov chains
  x <- mcmcchecks(x, autoburnin = autoburnin)
  
  # Define some notation
  Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
  Nvar <- ncol(x[[1]]) # number of variables
  Nchain <- length(x)

  if(Nvar == 1)stop("This function requires a Markov chain with multiple variables. If you have a univariate Markov chain, use asym.var() instead of asym.var.mat().")
  
  Ti <- lapply(x, getT, method = method, size = size)  # For each chain
  Tee <- matrix(Reduce("+", Ti)  / Nchain, nrow = Nvar)
  if(adjust == TRUE) Tee <- adjust_matrix(Tee, Niter)
  
  Tee

}

mcse.mat <- mcmcse:::mcse.mat
mcse.multi <- mcmcse:::mcse.multi


getT <- function(x, method, size) 
{
  mcse.multi(x, method = method, size = size)$cov
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


mcmcchecks <- function(x, autoburnin){
  
  # in case input is an mcmc object (single markov chain), change it to a matrix
  if(class(x) == "mcmc") {
    x <- as.matrix(x)
  }
  
  
  # in case input is a matrix (single markov chain), change it to a list
  if(class(x) == "matrix") {
    x <- list(x)
  }
  
  
  # in case input is of type mcmc.list, change it to a list of matrices
  if(class(x) == "mcmc.list") {
    x <- as.list(x)
    x <- lapply(x, as.matrix)
  }
  
  # make sure we have a list of matrices
  if(class(x) != "list") stop("Input x must be a list of matrices.")
  Nchain <- length(x) # number of chains. We also call this m.
  if(all.equal(lapply(x, class), as.list(rep("matrix", Nchain))) == FALSE) stop("Each item in list x must be a matrix.")
  
  # Define some notation.
  Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
  Nvar <- ncol(x[[1]]) # number of variables
  xnames <- colnames(x[[1]])
  
  # if multiple chains, ensure consistency in nrows, ncols
  if(Nchain > 1){    
    # check that each chain has some number of iterations
    if(all.equal(as.numeric(lapply(x, ncol)), rep(Nvar, Nchain)) != TRUE) stop("Unequal number of parameters between Markov chains. Each Markov chain must have the same number of columns.")
    
    # check that each chain has some number of iterations
    if(all.equal(as.numeric(lapply(x, nrow)), rep(Niter, Nchain)) != TRUE) stop("Unequal sample sizes between Markov chains. Each Markov chain must have the same number of rows.")
  }
  
  
  if (autoburnin && start(x) < end(x)/2) 
  {x <- window(x, start = end(x)/2 + 1)
  Niter <- nrow(x[[1]])
  }
  
  x
  
}

