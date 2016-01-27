`BiasCor` <-
function(b,D,ED,V,Vu,Vd){
# Input :
# b : the GLS estimated parameters
# D : the calculated design matrix
# ED : the expected value of the design matrix (when you calculate the design matrix in every place of your observed predictor variables you put in their ancestral value)
# V : the residual covariance matrix
# Vu : the predictor error covariance, can be either a matrix, vector of length(x) or a single number then we assume the matrix is diagonal with Vu on the diagonal 
# Vd : the predictor (phylogenetic) covariance matrix
# Output : the parameters bias corrected
    m<-length(b)    
    if (!is.matrix(Vu)){Vu<-diag(Vu,n,n)}
    V1<-solve(V)
    solve(diag(1,m,m)-pseudoinverse(t(D)%*%V1%*%D)%*%t(D)%*%V1%*%matrix(Vu%*%pseudoinverse(Vu+Vd)%*%(c(D)-c(ED)),ncol=ncol(D),nrow=nrow(D),byrow=FALSE))%*%b    
}

