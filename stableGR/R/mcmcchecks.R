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


