# x is a list of chains
# size is the batch mean size

size.and.trim <- function(x, size){
    
    # Define some notation
    Nchain <- length(x)
    
    if(is.vector(x[[1]])) #single component
    {
        Nvar <- 1
        Niter <- length(x[[1]])
    }
    if(is.matrix(x[[1]])) #univariate OR multivariate
    {
        Nvar <- ncol(x[[1]])
        Niter <- nrow(x[[1]])  # number of iterations per chains. We also call this n.
    } # number of variables
    
        
        
    if(is.null(size)){
        bvec <- sapply(x, FUN = batchSize, simplify = TRUE, method = "bm")  
        b <- floor(mean(bvec))
        a <- floor(Niter/b)
    }
    if(is.null(size) == FALSE){    # if size != NULL
        if (size == "sqroot") {   #calculation for square root
            b = floor(sqrt(Niter))
            a = floor(Niter/b)
        }  else if (size == "cuberoot") { #calculation for cube root
            b = floor(Niter^(1/3))
            a = floor(Niter/b)
        }  else {
            if (!is.numeric(size) || size <= 1 || size == Inf) 
                stop("'size' must be a finite numeric quantity larger than 1.")
            b = floor(size) #calculation for a user-specified batch size
            a = floor(Niter/b)
        }
    }
    
    ## trim away beginnings of each chain (if necessary)
    Nneeded <- a*b
    Ntrim <- Niter - Nneeded
    trim2 <- x
    if(Ntrim > 0){
        removethese <- 1:Ntrim
        # for(i in 1:Nchain){
        #   trimmedchain[[i]] <- trimchain(x[[i]], removethese)
        # }
        trim2 <- lapply(x, FUN = trimchain, removethese = removethese)
    }
    
    #return the size info and the list of trimmed chains
    list(a = a, b = b, Nneeded = Nneeded, Nvar = Nvar, trimmedchains = trim2) 
    
}

trimchain <- function(fullchain, removethese){
    Nvar <- ncol(fullchain)
    matrix(fullchain[-removethese,], ncol = Nvar)
}