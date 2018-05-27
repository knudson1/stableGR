\name{gelman.bm}
\alias{gelman.bm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gelman.bm(x, confidence = 0.95, transform = FALSE, autoburnin = TRUE, multivariate = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{confidence}{
%%     ~~Describe \code{confidence} here~~
}
  \item{transform}{
%%     ~~Describe \code{transform} here~~
}
  \item{autoburnin}{
%%     ~~Describe \code{autoburnin} here~~
}
  \item{multivariate}{
%%     ~~Describe \code{multivariate} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, confidence = 0.95, transform = FALSE, autoburnin = TRUE, 
    multivariate = TRUE) 
{
    x <- as.mcmc.list(x)
    if (nchain(x) < 2) 
        stop("You need at least two chains")
    if (autoburnin && start(x) < end(x)/2) 
        x <- window(x, start = end(x)/2 + 1)
    Niter <- niter(x)
    Nchain <- nchain(x)
    Nvar <- nvar(x)
    xnames <- varnames(x)
    if (transform) 
        x <- gelman.transform(x)
    x <- lapply(x, as.matrix)
    xbar <- matrix(sapply(x, apply, 2, mean, simplify = TRUE), 
        nrow = Nvar, ncol = Nchain)
    muhat <- apply(xbar, 1, mean)
    Si2 <- array(sapply(x, var, simplify = TRUE), dim = c(Nvar, 
        Nvar, Nchain))
    W <- apply(Si2, c(1, 2), mean)
    Ssq <- diag(W)
    s2 <- matrix(apply(Si2, 3, diag), nrow = Nvar, ncol = Nchain)
    var.s2 <- apply(s2, 1, var)/Nchain
    tau2i <- sapply(x, gettau)
    tau2 <- apply(tau2i, 1, mean)
    tau2var <- apply(tau2i, 1, var)/Nchain
    cov.s2t2 <- sapply(1:Nvar, function(index) {
        cov(x = tau2i[index, ], y = s2[index, ])
    })/Nchain
    sigsq <- (Niter - 1) * Ssq/Niter + tau2/Niter
    V <- sigsq + tau2/(Niter * Nchain)
    var.V <- ((Niter - 1)/Niter)^2 * var.s2/Nchain + +((Nchain + 
        1)/(Nchain * Niter))^2 * tau2var/Nchain + +cov.s2t2/Nchain
    df.V <- (2 * V^2)/var.V
    df.adj <- (df.V + 3)/(df.V + 1)
    B.df <- Nchain - 1
    W.df <- (2 * Ssq^2)/var.s2
    R2.fixed <- (Niter - 1)/Niter
    R2.random <- (1 + 1/Nchain) * (1/Niter) * (tau2/Ssq)
    R2.estimate <- R2.fixed + R2.random
    R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * 
        R2.random
    psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
    dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
    out <- list(psrf = psrf)
    class(out) <- "gelman.diag"
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
