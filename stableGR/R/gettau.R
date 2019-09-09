gettau <- function(x1, method, size) 
{
	(mcse.mat(x1, method = method, size = size)[ ,2])^2 
}

