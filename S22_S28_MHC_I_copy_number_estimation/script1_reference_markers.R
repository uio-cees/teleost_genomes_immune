###################################################################
###--- Script 1 of the procedure.
###---
###--- Here we first decide which reference markers are possible
###--- to use for estimating copy numbers.
###--- Step 1 - Sanity check on read coverage.
###--- Step 2 - Check on effect of BLAST-cutoff.
###--- Step 3 - Final check if copy number estimates are not 0.
###################################################################
options( stringsAsFactors = FALSE )
source( "copynumfun.R" )



### Loads data and print som summary
cat( "Loading data...\n" )
load( "data/Teleost_all_2015.RData" )
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



### First data sanity check - checking for too low coverage
### Some reference markers may, for some reason, have extremely low coverage
### We would like to exclude these since they may affect all other estimates
### of copy numbers later
low.cov <- matrix( 0, nrow=n.fish, ncol=n.ref )
rownames( low.cov ) <- colnames( ref.hits )
colnames( low.cov ) <- uref
coverage <- (N.reads*marker.length)/G.size    # This is the expected coverage per fish
for( i in 1:n.ref ){
  cat( "Reference ", uref[i], "\n" )
  
  idd <- seq( 1, nrow(ref.hits), 4 )                          # the hits under cutoff -5
  phi1 <- mean( ref.hits[idd[i],]/coverage )                  # very rough estimate of cutoff-factor
  ehits <- coverage*phi1                                      # expected number of hits based on phi1
  rr1 <- sign( ref.hits[idd[i],] - (ehits-3*sqrt(ehits)) )    # if rr1 is negative the coverage is very low
  
  idd <- seq( 2, nrow(ref.hits), 4 )                          # the hits under cutoff -10
  phi2 <- mean( ref.hits[idd[i],]/coverage )                  # very rough estimate of cutoff-factor
  ehits <- coverage*phi2                                      # expected number of hits based on phi2
  rr2 <- sign( ref.hits[idd[i],] - (ehits-3*sqrt(ehits)) )    # if rr2 is negative the coverage is very low

  idd <- seq( 3, nrow(ref.hits), 4 )                          # etc...
  phi3 <- mean( ref.hits[idd[i],]/coverage )
  ehits <- coverage*phi3
  rr3 <- sign( ref.hits[idd[i],] - (ehits-3*sqrt(ehits)) )

  idd <- seq( 4, nrow(ref.hits), 4 )
  phi4 <- mean( ref.hits[idd[i],]/coverage )
  ehits <- coverage*phi4
  rr4 <- sign( ref.hits[idd[i],] - (ehits-3*sqrt(ehits)) )
  
  rmat <- matrix( c(rr1,rr2,rr3,rr4), ncol=4, byrow=F )       # The rr1,...,rr4 values
  low.cov[,i] <- rowSums( rmat )                              # summing the signs, i.e. need 3 out of 4 negative to
}                                                             # get a negative sum
ref.keep <- (low.cov >= 0 )
f.disc <- 
r.disc <- colSums( !ref.keep )
cat( "Number of discarded reference markers per fish:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",rowSums(!ref.keep)[i],"reference markers\n")
} )
cat( "Number of discarded fishes per reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",colSums( !ref.keep )[i],"fishes\n")
} )




### Manual curation
### These two references were discarded in the end because they are
### probably too closely linked (physical position) to some other references
ref.keep[,2] <- FALSE  # Ref 1b is OUT
ref.keep[,7] <- FALSE  # Ref 2d is OUT




### Next, checking which reference genes have read-hits that
### cannot be explained by the cutoff-model.
### In short, we expect that as the blast-cutoff is made stricter, the
### number of hits should decrease in a smooth way.
### Here we look for cases where this is clearly violated
mse <- matrix( 0, nrow=n.fish, ncol=n.ref )
colnames( mse ) <- uref
rownames( mse ) <- colnames( ref.hits )
for( ss in 1:n.fish ){
  cat( rownames(mse)[ss], ":\n" )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  idx.in <- which( ref.keep[ss,] )
  n.in <- length( idx.in )
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  y.mat <- y.mat[,idx.in]
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  
  g.mat <- matrix( rep( lst$Gamma, 4 ), nrow=4, byrow=TRUE )
  p.mat <- matrix( rep( lst$Phi, n.in ), ncol=n.in, byrow=FALSE )
  y.hat <- C * g.mat * p.mat
  r.mat <- y.mat - y.hat
  mse[ss,idx.in] <- apply( r.mat, 2, function(x){sum(x^2)/(4-1)} )
  rr <- range( c( as.vector( y.mat ), as.vector(y.hat) ) )
  plot( rr, rr, type="l", col="red", main=rownames(mse)[ss], 
        xlab="Observed number of reads", ylab="Predicted number of reads" )
  points( y.mat, y.hat, pch=16 )
  Sys.sleep( 0.2 )
}
sigma2 <- mean( mse, trim=0.01 )
lambda <- (4-1)*mse/sigma2
limit <- qchisq( 0.99, df=(4-1) )
ref.keep <- (lambda <= limit)&ref.keep
cat( "Number of discarded reference markers per fish:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",rowSums(!ref.keep)[i],"reference markers\n")
} )
cat( "Number of discarded fishes per reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",colSums( !ref.keep )[i],"fishes\n")
} )




### Manual curation again
### Fish 21 suffers from having few reference markers in the 
### analysis above. Inspect the plot for fish_21, it looks fine!
### We decide to keep reference 8 and 13 even if they are
### excluded in the procedure above.
ref.keep[21,14] <- TRUE  # Ref 8 is IN
ref.keep[21,20] <- TRUE  # Ref 13 is IN





### Third, eliminating reference genes that still estimates 
### to 0 copies, as these will blow the other copy number 
### estimates sky high later...
ref.hat <- matrix( 0, nrow= n.fish, ncol=n.ref )
for( ss in 1:n.fish ){
  cat( rownames(ref.keep)[ss], ":\n" )
  C <- (N.reads[ss]*marker.length)/G.size[ss]
  idx.in <- which( ref.keep[ss,] )
  n.in <- length( idx.in )
  y.mat <- matrix( ref.hits[,ss], nrow=4, ncol=n.ref, byrow=F )
  y.mat <- y.mat[,idx.in]
  lst <- copynum.iter.iter( y.mat, C, rep( 1, n.in ), as.integer=FALSE )
  ref.hat[ss,idx.in] <- lst$Gamma
}
ref.keep <- (round(ref.hat) > 0)&ref.keep
cat( "Number of discarded reference markers per fish:\n" )
x <- sapply( 1:n.fish, function(i){
  cat(rownames(ref.keep)[i], "discards",rowSums(!ref.keep)[i],"reference markers\n")
} )
cat( "Number of discarded fishes per reference marker:\n" )
x <- sapply( 1:n.ref, function(i){
  cat(colnames(ref.keep)[i], "is discarded in",colSums( !ref.keep )[i],"fishes\n")
} )


save( ref.keep, file="res/ref_keep.RData" )
