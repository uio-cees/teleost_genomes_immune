`ResidV` <-
function(beta, topology, times, alpha, vy, response, me.response, random.cov, me.random.cov, mecov.random.cov=NULL, ultrametric=TRUE)
{
#### will need to do some precalculations here
a<-alpha
Y <- response[!is.na(response)]; 
N <- length(Y);
pred<-data.frame(random.cov);
n.pred<-length(pred[1,]);
pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
if(is.null(me.response)) me.response<-diag(rep(0, times=length(response[!is.na(response)])))  else me.response<-diag(me.response[!is.na(me.response)]);
if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);
T <- times[terminal.twigs(topology)];
tia<-tsia(topology, times);
tja<-tsja(topology, times);
term<-terminal.twigs(topology);
pt<-parse.tree(topology, times);
ta<-pt$bt;
tij<-pt$dm;
num.prob<-matrix(data=0, nrow=N, ncol=N) #this matrix is included for cases where species split at the root;
cm2<-matrix(data=0, nrow=N, ncol=N);
  s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
  for(i in 1:n.pred)
   {
    s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2]);
   }
 beta1<-beta
if(ultrametric==TRUE)
 s1<-as.numeric(s.X%*%((beta1[2:(n.pred+1),])*(beta1[2:(n.pred+1),]))) 

else
 s1<-as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))   


for(p in 1:N)
 {
  for(q in 1:N)
   {
     if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p])
   }
 }
 
cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij)
 
 for(p in 1:N)
  {
   for(q in 1:N)
    {
     cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]))
    }
  }
  
 if(ultrametric==TRUE)
  {
   mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[2:(n.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
   mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))  
   V<-cm1+(s1*ta*cm2)+me.response+mv-mcov 
  }
 else
  {
   mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[4:(n.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))  
   mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))   
   V<-cm1+(s1*ta*cm2)+me.response+mv-mcov      
  }
  return(V)  
}

