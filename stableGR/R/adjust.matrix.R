adjust.matrix <- function(mat, N, epsilon = sqrt(log(N)/dim(mat)[2]), b = 1/2)
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

