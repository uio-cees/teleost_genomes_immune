`np.regression` <-
function(response, me.response=NULL, random.cov, me.random.cov=NULL, mecov.random.cov=NULL, convergence=NULL){
  if(is.null(convergence)) convergence=0.000001 
  Y <- response[!is.na(response)];
  N <- length(Y);
	n.fixed=1;
  pred<-data.frame(random.cov)
  n.pred<-length(pred[1,])
  pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred)	
 
	
## Code for optional measurement variances	
if(is.null(me.response)) me1<-diag(rep(0, times=length(response[!is.na(response)])))  else me1<-diag(me.response[!is.na(me.response)]);
if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);

  
  
  X<-cbind(1, pred)
  x.ols<-cbind(1, pred)
  V1<-diag(rep(1, times=N))
  beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
  eY<-X%*%beta1;  
  pred.mean<-X%*%beta1
     g.mean<-(t(rep(1, times=N))%*%solve(V1)%*%Y)/sum(solve(V1));
     sst<-t(Y-g.mean)%*% solve(V1)%*%(Y-g.mean)
     sse<-t(Y-pred.mean)%*%solve(V1)%*%(Y-pred.mean)
     sigma<-sse/(N-(n.pred+1))
   repeat{
    V<-diag(rep(sigma, times=N))+me1 + diag(as.numeric(me.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])));
    V.inverse<-solve(V)
    beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
    test<-matrix(nrow=(n.pred+1))
     for(f in 1:(n.pred+1))
        {
         if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break  
      beta1<-beta.i
   }
   beta1<-beta.i
    beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
   eY<-X%*%beta1;  
   pred.mean<-X%*%beta1
     g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
     sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
     sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
     sigma<-sse/(N-(n.pred+1))
     r.squared<-(sst-sse)/sst 

   sigma.sq.y<-(1/N)*(sum((Y-X%*%beta1)^2))
   log.like<-(-N/2)*log(2*pi)-(N/2)*log(sigma.sq.y)-sum((Y-X%*%beta1)^2/(2*sigma.sq.y))
   ml<-log.like
   reg.par<-matrix(data=0, nrow=(n.pred+1), ncol=2, dimnames=list(c("Intercept", if(n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
   reg.par[,1]<-round(beta1,5);
   reg.par[,2]<-round(sqrt(diag(beta.i.var)),5)
   modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
   n.par=n.pred;
   modfit[1,1]=ml 
   modfit[2,1]=-2*ml+2*(2+n.par)
   modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
   modfit[4,1]=-2*ml+log(N)*(2+n.par)
   modfit[5,1]=r.squared*100 
   modfit[6,1]=sst
   modfit[7,1]=sse
   message("REGRESSION PARAMETERS");message("");
   print(reg.par);message("");
   message("--------------------------------------------------"); 
   message("MODEL FIT");message("");
   print(modfit); message("");
   message("==================================================");
}

