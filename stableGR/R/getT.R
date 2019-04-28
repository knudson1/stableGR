getT <- function(x, method, size) 
{
  mcse.multi(x, method = method, size = size)$cov
}

