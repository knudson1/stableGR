gettau <- function(x1, method) 
{
	(mcse.mat(x1, method = method)[ ,2])^2 
}

