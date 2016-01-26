


copynum.iter.iter <- function( y.mat, C, prior, max.iter=10, as.integer=TRUE, verbose=TRUE ){
  # y.mat must have one row for each stringency and one column for each gene
  gamma.hat <- prior
  gamma.old <- rep( 1000, length( prior ) )
  if( verbose ) cat( "copynum.iter.iter:\n" )
  tol <- 0.01
  
  niter <- 0
  while( (max( abs( gamma.hat - gamma.old ) ) > tol) & (niter<max.iter) ){
    gamma.old <- gamma.hat
    lst <- copynum.iter( y.mat, C, prior, max.iter, as.integer, verbose )
    gamma.hat <- lst$Gamma
    phi.hat <- lst$Phi
    prior <- round( gamma.hat )
    niter <- niter + 1
    if( verbose ) cat( " iteration.iteration", niter, "\n" )
  }
  return( list( Gamma=gamma.hat, Phi=phi.hat ) )
}


copynum.iter <- function( y.mat, C, prior, max.iter=10, as.integer=TRUE, verbose=TRUE ){
# y.mat must have one row for each stringency and one column for each gene
  dd <- dim( y.mat )
  phi.hat <- rep( 1, dd[1] )
  gamma.hat <- prior
  gamma.old <- rep( 1000, dd[2] )
  if( verbose ) cat( "copynum.iter:\n" )
  tol <- 0.01
  
  niter <- 0
  while( (max( abs( gamma.hat - gamma.old ) ) > tol) & (niter<max.iter) ){
    gamma.old <- gamma.hat
    g.mat <- matrix( rep( gamma.hat, dd[1] ), nrow=dd[1], byrow=T )
    idx <- which( g.mat <= 0.5, arr.ind=TRUE )
    g.mat[idx] <- NA
    phi.hat <- rowMeans( y.mat/g.mat, na.rm=TRUE )/C
    p.mat <- matrix( rep( phi.hat, dd[2] ), nrow=dd[1], byrow=F )
    gamma.hat <- colMeans( y.mat/p.mat )/C
    if( as.integer ) gamma.hat <- round( gamma.hat )
    niter <- niter + 1
    if( verbose ) cat( " iteration", niter, "\n" )
  }
  return( list( Gamma=gamma.hat, Phi=phi.hat ) )
}


copynum <- function( y, ref.mat, ref.gamma, as.integer=TRUE ){
  dd <- dim( ref.mat )
  idx <- which( ref.mat == 0, arr.ind=T )
  ref.mat[idx] <- 1
  y.mat <- matrix( rep( y, dd[2] ), nrow=dd[1], byrow=F )
  g.mat <- matrix( rep( ref.gamma, dd[1] ), nrow=dd[1], byrow=T )
  gamma.hat <- g.mat * y.mat/ref.mat
  if( as.integer ) gamma.hat <- round( gamma.hat )
  gamma.hat[which(gamma.hat==Inf,arr.ind=T)] <- NA
  return( gamma.hat )
}