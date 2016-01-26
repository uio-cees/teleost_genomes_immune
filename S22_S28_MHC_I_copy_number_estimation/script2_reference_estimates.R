###################################################################
###--- Script 2 of the procedure.
###---
###--- You need to run script1 prior to this.
###--- Here we estimate the reference markers copy numbers
###--- Step 1 - Point estimates
###--- Step 2 - Bootstrapping 
###################################################################
options( stringsAsFactors = FALSE )
source( "copynumfun.R" )




### Loads data and print som summary
cat( "Loading data...\n" )
load( "data/Teleost_all_2015.RData" )
load( "res/ref_keep.RData" )
uref <- unique( gene )
uref <- uref[3:length(uref)]             # two first are not reference markers
r.rows <- 9:nrow( hits )                 # eight first row are not reference markers
ref.hits <- hits[r.rows,]
marker.length <- 90*3                    # marker length in bases
G.size <- colMeans( genome.size )        # eight different estimates of genome size, using the mean
n.fish <- ncol( ref.hits )
n.ref <- length( uref )
cuts <- c(-5,-10,-15,-20)
cat( "   we have data for", n.ref, "reference markers\n" )
cat( "   we have data for", n.fish, "fishes\n" )
cat( "   reference markers are", marker.length, "bases long\n" )
x <- sapply( 1:n.fish, function(i){
  cat(colnames(ref.hits)[i], "has", G.size[i],"basepairs and",N.reads[i],"reads\n")
} )





### Point estimates of reference marker copy number in each fish
### based on read-counts and ALL reference markers for the fish.
### The matrix ref.point contains one row for each fish (66 rows)
### and one column for each marker (21), but since some
### marker have been excluded for the various fish (see script1)
### there are some NA in the matrix.
ref.point <- matrix( NA, nrow=ncol(ref.hits), ncol=n.ref )
rownames( ref.point ) <- names( G.size )
colnames( ref.point ) <- uref
for( ss in 1:n.fish ){
  cat( rownames(ref.point)[ss], ":\n" )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  idx.in <- which( ref.keep[ss,] )
  n.in <- length( idx.in )
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  y.mat <- y.mat[,idx.in]
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  ref.point[ss,idx.in] <- lst$Gamma
}
save( ref.point, file="res/ref_point.RData" )
boxplot(ref.point,las=2,pch=16,ylab="Reference marker copy number")






### The bootstrap procedure
### The difference to the point estimate above is that
### estimates are based on a bootstrap-sample of the reference
### markers (not all markers), and this is repeated 1000
### times.
### When estimating the copy number for marker g, this marker must of course
### be included, but the OTHER markers are bootstrapped 1000 times, and this
### is repeated for every marker for each fish.
### The 3-dimensional array ref.boot stores the results.
phi.hat <- matrix( -1, nrow=4, ncol=n.fish )
rownames( phi.hat ) <- paste( "log10(E)=", cuts, sep="" )
colnames( phi.hat ) <- colnames( ref.hits )
N.boot <- 1000
ref.boot <- array( NA, dim=c(n.fish, n.ref, N.boot), dimnames=list( Fish=names(G.size), REF=uref, Boot=1:N.boot ) )
for( ss in 1:n.fish ){
  cat( rownames(ref.point)[ss] )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  idx.in <- which( ref.keep[ss,] )
  y.mat <- y.mat[,idx.in]
  n.in <- length( idx.in )
  for( g in 1:n.in ){
    idb <- which( 1:n.in != g )
    for( b in 1:N.boot ){
      idx.boot <- c( g, sample( idb, size=(n.in-1), replace=TRUE ) )
      y.boot <- y.mat[,idx.boot]
      lst <- copynum.iter.iter( y.boot, C, rep( 1, n.in ), max.iter=10, as.integer=FALSE, verbose=FALSE )
      ref.boot[ss,idx.in[g],b] <- lst$Gamma[1]
    }
    cat( "." )
  }
  cat( "\n" )
}

save( ref.boot, file="res/ref_boot.RData" )

