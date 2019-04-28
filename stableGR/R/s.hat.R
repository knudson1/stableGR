s.hat <- function(x){

  # Define some notation.
  Nvar <- ncol(x[[1]]) # number of variables
  Nchain <- length(x)
  
	Si2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, Nvar, Nchain)) 
	W <- apply(Si2, c(1, 2), mean) # Average the vcov matrices across chains. #aka S
	W
}
