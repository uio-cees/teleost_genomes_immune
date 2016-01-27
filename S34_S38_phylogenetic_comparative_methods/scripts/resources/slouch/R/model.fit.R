`model.fit` <-
function(topology, times, half_life_values, vy_values, response, me.response=NULL, fixed.fact=NULL,fixed.cov=NULL, me.fixed.cov=NULL, mecov.fixed.cov=NULL, random.cov=NULL, me.random.cov=NULL, mecov.random.cov=NULL,  intercept="root", ultrametric=TRUE, support=NULL, convergence=NULL)
{
  
# SET DEFAULTS IF NOT SPECIFIED	
	
 if(is.null(support)) support=2;
 if(is.null(convergence)) convergence=0.000001;
 if(is.null(me.response)) me.response<-diag(rep(0, times=length(response[!is.na(response)])))  else me.response<-diag(me.response[!is.na(me.response)]);

 
# DETERMINE MODEL STRUCTURE FROM INPUT AND WRITE A SUMMARY TO THE R CONSOLE
 
 if(is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "IntcptReg";
 if(!is.null(fixed.fact) && is.null(fixed.cov) && is.null(random.cov)) model.type <- "ffANOVA";
 if(!is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <-"ffANCOVA";
 if(!is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmANCOVA";
 if(!is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mmfANCOVA";   if(is.null(fixed.fact) && is.null(fixed.cov) && !is.null(random.cov)) model.type <- "rReg";
 if(is.null(fixed.fact) && !is.null(fixed.cov) && is.null(random.cov)) model.type <- "fReg";
 if(is.null(fixed.fact) && !is.null(fixed.cov) && !is.null(random.cov)) model.type <- "mfReg"; 
 
# Write type of model to screen
 
 message("")
 message("MODEL SUMMARY")
 message("")
 if(model.type=="IntcptReg")
 {
 	 message("You have specified an OU model for a response variable regressed on a grand mean, i.e. one global optima");
 	 if(ultrametric==FALSE) 
 	 {
 	 GS_head<-c("Ya", "Theta_Global"); 
 	 n.par<-2;
 	 }
 	 else
 	 { 
 	 GS_head<-("Theta_Global");
 	 n.par<-1;
 	 }
 }
 else
 if(model.type=="ffANOVA" ) 
  {
     message("You have specified an OU model for a response variable modeled on optima determined by fixed, categorical predictor variables");
     if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
  }
 
 else
 if(model.type=="ffANCOVA")
  {
   message("You have specified an OU model for a response variable modeled on optima determined by both fixed categorical predictors and an instantaneous scaling with a fixed covariate"); 
   if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));
   
   }
 
 
 else

 if(model.type=="mmANCOVA") 
  {
   message("You have specified an OU model for a response variable modeled on optima determined by both fixed, categorical factors as well as covariates which themselves randomly evolve (modeled as Brownian-motions)");
   
   
   if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));  
  }
  
   if(model.type=="mmfANCOVA") 
  {
   message("You have specified an OU model for a response variable modeled on optima determined by both fixed, categorical factors as well as covariates which themselves randomly evolve (modeled as Brownian-motions)");
   
   
   if(is.null(intercept)) GS_head<-c("Ya", levels(as.factor(fixed.fact))) else GS_head<-levels(as.factor(fixed.fact));  
  }
  else
  
  

 if(model.type=="rReg") message("You have specified an OU model for a response variable modeled on optima that are determined by randomly evolving covariates (modeled as Brownian-motions)")
 
 else
 if(model.type=="fReg") message("You have specified an OU model for a response variable modeled on optima that are determined by an instantaneous scaling with fixed covariates")
 
  else
 if(model.type=="mfReg") message("You have specified an OU model for a response variable modeled on optima that are determined by both an instantaneous scaling with fixed covariates and randomly evolving covariates (modeled as Brownian-motions)"); 
 message("")
 
# Summarize dataset, response, predictors,  tree height and sample size and write to screen
 
ms<-list(Dataset=search()[2], Response=deparse(substitute(response)), Fixed.factor=deparse(substitute(fixed.fact)),Fixed.covariates=deparse(substitute(fixed.cov)), Random.covariates=deparse(substitute(random.cov)), Sample.size=length(response[!is.na(response)]), Tree.height=max(times))
ms<-as.matrix(ms)
colnames(ms)<-"Summary"
print(ms)
message("")
message("GRID SEARCH PARAMETER SUPPORT")
message("")

# SPECIFY COMPONENTS THAT ARE COMMON TO ALL MODELS

Y <- response[!is.na(response)]; 
N <- length(Y);
T <- times[terminal.twigs(topology)];
tia<-tsia(topology, times);
tja<-tsja(topology, times);
term<-terminal.twigs(topology);
pt<-parse.tree(topology, times);
ta<-pt$bt;
tij<-pt$dm;
num.prob<-matrix(data=0, nrow=N, ncol=N) #this matrix is included for cases where species split at the root;
cm2<-matrix(data=0, nrow=N, ncol=N);
gof<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values), dimnames=list(half_life_values, vy_values));
h.lives<-matrix(data=0, nrow=length(half_life_values), ncol=length(vy_values))
ln2<-log(2)
half_life_values<-rev(half_life_values)



# EVALUATE IF IT IS A FIXED FACTOR PREDICTOR OR INTERCEPT ONLY MODEL THEN SET UP APPROPRIATE DESIGN AND VARIANCE MATRICES AND ESTIMATE PARAMETERS WITHOUT ITERATED GLS

if(model.type =="IntcptReg" || model.type == "ffANOVA") 
 {
 if(model.type=="IntcptReg") regime.specs<-rep(1, times=length(topology)) else regime.specs<-fixed.fact;

 cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head), sep="   ");  
 message(" ");
 
  
  for(i in 1:length(half_life_values))
   {
    for(k in 1:length(vy_values))
     { 
      vy <- vy_values[k];
      if(half_life_values[i]==0)
       {
        a<-1000000000000000000000; 
        V<-diag(rep(vy, times=N)) + me.response;
       }
      else 
       {
       	a <- ln2/half_life_values[i];
       	V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response;
       }
       
       if(model.type=="IntcptReg")
       {
       	if(half_life_values[i]==0 ||a>=1000000000000000000000) X<-matrix(data=1, nrow=N, ncol=1)       	    else   	
       if(ultrametric==TRUE) X<-matrix(data=1, nrow=N, ncol=1) 
       else
        {
         X<-matrix(data=0, nrow=N, ncol=2);
         X[,1]<-1-exp(-a*T);
         X[,2]<-exp(-a*T)
       	}
       }
       else       
       X<-weight.matrix(a, topology,times, N, regime.specs, fixed.cov, intercept) 
                 
    # GLS estimation of parameters for fixed model
    
 V.inverse<-solve(V)
 beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
 beta0<-beta.i
 eY<-X%*%beta0
 resid<-Y-eY				
 gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
 print(c(half_life_values[i], vy, round(gof[i,k], 4), round(as.numeric(t(beta0)), 4)))

       
     } # end of half-life loop
    } # end of vy loop
    
# Search GOF matrix for best estimates of alpha and vy #

   x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
    for(i in 1:length(half_life_values))
     {
      for(j in 1:length(vy_values))
       {
        if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
       }
     }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support) gof[i, j]=ml-support;
      }
    }     

   gof=gof-ml
   
 # final GLS estimations for corrected optima using best alpha and vy estimates #

  if(alpha.est==Inf) alpha.est<-1000000000000000000000 
  
  if(model.type=="IntcptReg")
       {
       	if(alpha.est==Inf || alpha.est>=1000000000000000000000 ) X<-matrix(data=1, nrow=N, ncol=1)  
       	else
       	       if(ultrametric==TRUE) X<-matrix(data=1, nrow=N, ncol=1) 
       else
        {
         X<-matrix(data=0, nrow=N, ncol=2);
         X[,1]<-1-exp(-alpha.est*T);
         X[,2]<-exp(-alpha.est*T)
       	}
       }
       else       
     X<-weight.matrix(alpha.est, topology,times, N, regime.specs, fixed.cov, intercept)  
  
  V<-((vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)) + me.response;
  V.inverse<-solve(V);
  beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X);
  beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y);
  gls.beta0<-beta.i; 
  
 # code for calculating SSE, SST and r squared
 
  pred.mean <- X%*%gls.beta0
  g.mean <- (t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst <- t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse <-t (Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared <- (sst-sse)/sst  



      
 } # END OF FIXED PREDICTOR OR INTERCEPT ONLY PARAMETER ESTIMATION
 
 
if(model.type =="ffANCOVA" || model.type == "fReg")  
{
	
 fixed.pred<-data.frame(fixed.cov);
 n.fixed.pred<-length(fixed.pred[1,]);
 fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
 if(is.null(me.fixed.cov)) me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred);
 if(is.null(mecov.fixed.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.cov<-matrix(data=me.cov.fixed.cov[!is.na(me.cov.fixed.cov)], ncol=n.fixed.pred);  
 
if(model.type=="fReg")
  {
  	x.ols<-cbind(1, fixed.pred);
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y); 
    n.fixed<-1
    cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Bo", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), sep="   ");  
    message("");	      
  }
  
if(model.type=="ffANCOVA") 
  {
  	regime.specs<-fixed.fact;
  	n.fixed<-length(levels(as.factor(regime.specs)))
    regime.specs<-as.factor(regime.specs)
    x.ols<-weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept);  	beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);
    cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), sep="   "); 
    message("");  
   }
   
   
for(i in 1:length(half_life_values))
  {  
   for(k in 1:length(vy_values))
    { 
      vy <- vy_values[k];	
      if(half_life_values[i]==0) 
       {
        a<-1000000000000000000000
        V<-diag(rep(vy, times=N)) + me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 
       
       }
      else 
       {
       	a <- ln2/half_life_values[i];
       	V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 
       }
      if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept); 
      
      ##### iterated GLS 
      
       con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
        if(half_life_values[i]==0)
         {
          a<-1000000000000000000000           	       
          V<-diag(rep(vy, times=N)) + me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 
                    
         }
        else
         {	
          a <- ln2/half_life_values[i];
       	  V<-((vy)*(1-exp(-2*a*ta))*exp(-a*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 

         } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT
         
         if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept); 
         
         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #   
         
         V.inverse<-solve(V)
                   
          beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
         }  
eY<-X%*%beta1
 resid<-Y-eY
 gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
 print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4))) 

 
 ### END OF ITERATED GLS
       
     } # end of half-life loop
    } # end of vy loop
    
 # Search GOF matrix for best estimates of alpha and vy #

   x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
    for(i in 1:length(half_life_values))
     {
      for(j in 1:length(vy_values))
       {
        if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
       }
     }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support) gof[i, j]=ml-support;
      }
    }     

   gof=gof-ml
     
   # final GLS estimations for corrected optima using best alpha and vy estimates #

 con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
        if(alpha.est==Inf)
         {
          a<-1000000000000000000000           	       
          V<-diag(rep(vy, times=N)) + me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 
                    
         }
        else
         {	
          
       	  V<-((vy)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij))+me.response + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 

         } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT
         
         if(model.type=="fReg") X<-cbind(1, fixed.pred) else  X<-weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept); 
         
         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #   
         
         V.inverse<-solve(V)
                   
          beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
         }  
         
       gls.beta0<-beta1;
       beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X);
           
 # code for calculating SSE, SST and r squared
 
  pred.mean <- X%*%gls.beta0
  g.mean <- (t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst <- t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse <-t (Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared <- (sst-sse)/sst  
      
 	
 }  # END OF fReg AND ffANCOVA ESTIMATION ROUTINES must still add iterated GLS for me
 
 
 
 # EVALUATE IF IT IS A FIXED MODEL ANCOVA, MIXED MODEL ANCOVA OR RANDOM PREDICTOR REGRESSION, ESTIMATE PARAMETERS WITH ITERATED GLS TO A) TAKE MEASUREMENT VARIANCE INTO ACCOUNT OR B) RANDOM EFFECTS INTO ACCOUNT IN THE CASE OF THE MIXED MODEL AND REGRESSION
    
if(model.type == "mmANCOVA" || model.type=="rReg")  ### more models here
{
 # SET UP INITIAL MATRICES FOR MULTIPLE REGRESSION AND CALCULATE THETA AND SIGMA FOR RANDOM PREDICTOR / S
 
 pred<-data.frame(random.cov);
 n.pred<-length(pred[1,]);
 pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
  if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
 if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);
  
 s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
   for(i in 1:n.pred)
    {
    s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2]);
    }
      
   theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
   for(i in 1:n.pred)
    {
     theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1]);
    }    
 
 # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
 
  
## INITIAL OLS ESTIMATES TO SEED ITERATED GLS

 if(model.type=="rReg")  
  {
  	x.ols<-cbind(1, pred);
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);  
    if(ultrametric == FALSE) beta1<-rbind(0, 0, beta1); # 2 additional parameter seeds for Ya and Xa    
  }
  
  if(model.type=="mmANCOVA")
  {
  	regime.specs<-fixed.fact;
  	n.fixed<-length(levels(as.factor(regime.specs)))
    regime.specs<-as.factor(regime.specs)
    x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred);  	beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);  
  }	
  
 
  
  # GRID ESTIMATION ROUTINE AND ITERATED GLS FOR MODELS THAT INCLUDE RANDOM EFFECTS
  
 if(model.type=="mmANCOVA")   
 {
  
 cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   "); 
 
 
 message(" ");	
 	
 for(i in 1:length(half_life_values))
  {  
   for(k in 1:length(vy_values))
    { 
      if(half_life_values[i]==0) a<-1000000000000000000000 else a <- ln2/half_life_values[i];
      vy <- vy_values[k];
      X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-a*T))/(a*T))*pred); 
       if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1)); n.fixed<-n.fixed+1}
       if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}        
                                  
    # CODE FOR ESTIMATING BETA USING ITERATED GLS 
    
       con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
        if(half_life_values[i]==0)
          {
           X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred);    
           V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),]))); 
          }
        else
         {	
         	
          X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
          
                                              	
          s1<-as.numeric(s.X%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),]));   
          
          
          for(p in 1:N)
           {
            for(q in 1:N)
             {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
             }
           } 
          cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
          for(p in 1:N)
           {
            for(q in 1:N)
             {
             cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/ (a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
             }
           }

      mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));   
           mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));  
           
                      
           V<-cm1+(s1*ta*cm2)+me.response+mv-mcov  
          } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #   
         
         V.inverse<-solve(V)
         if(half_life_values[i]==0)
          {
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
         else
          {
         
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
          }
         }  
                  
             
      ### END OF ITERATED GLS ESTIMATION FOR BETA #
      
     if(half_life_values[i]==0)
      {
       X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)   
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1):length(beta1),])))  
       V.inverse<-solve(V)
       eY<-X%*%beta1
       resid<-Y-eY;
       gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }
      else
      {
       s1<-as.numeric(s.X%*%(beta1[(n.fixed+1):length(beta1),]*beta1[(n.fixed+1):length(beta1),]))  
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
         }
       } 
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
         }
       }
           X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-a*T))/(a*T))*pred);          
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));  
         V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
        V.inverse<-solve(V)
     
       eY<-X%*%beta1
      
      resid<-Y-eY; 
      gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }  # END OF CONDITION FOR HALF-LIFE = 0 #
      print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4))) 
     }
    }
    
    
    
# END OF GRID SETUP,START OF GRID SEARCH FOR BEST ALPHA AND VY ESTIMATES #

   x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
      }
    }
   gof=gof-ml
   
  
   n.fixed<-length(levels(as.factor(regime.specs)))   ### reset before final regression
   
   
  # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #
   
  if(alpha.est==Inf || alpha.est >=1000000000000000000000)
   {
   	x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)   
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y) 
    con.count<-0;
    repeat
   	 {
      s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),]))
      X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)  
      V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),]))) 
      V.inverse<-solve(V)
      beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(length(beta.i)))
       for(f in 1:(length(beta.i)))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break 
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
         gls.beta1<-beta.i
       }
     gls.beta1<-beta.i
     X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)  
     V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1):length(beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1):length(gls.beta1),]))) 
     pred.mean<-X%*%gls.beta1
     g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
     sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
   
     sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    
     r.squared<-(sst-sse)/sst 
         }
   else
   {
    x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)   
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y) 
    con.count<-0;
    
    X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred); 
       if(length(X[1,]) > length(gls.beta1)) {gls.beta1<-as.matrix(c(0, gls.beta1)); n.fixed<-n.fixed+1}
       if(length(X[1,])< length(gls.beta1)) {gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs)))}        
   repeat
   {
   	
   	X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
    s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1):length(gls.beta1),]*gls.beta1[(n.fixed+1):length(gls.beta1),]))  
   	     
 
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
    
                 
            mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))    
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))*2), ncol=n.pred)))    
       V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
          V.inverse<-solve(V)
           beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(length(beta.i)))
       for(f in 1:(length(beta.i)))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break  
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

          
         gls.beta1<-beta.i
         
          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.cov, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
            
                mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))   
     mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(gls.beta1[(n.fixed+1):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))*2), ncol=n.pred)))   
    V<-cm1+(s1*ta*cm2)+me.response+mv-mcov; 
    pred.mean<-X%*%gls.beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst   
       }  
      }     

                # END OF ITERATED GLS LOOP #


            
                
      	
 } # END OF ESTIMATION MIXED MODEL ANCOVA	
 
 
if(model.type=="rReg")   
 {
 
 if(ultrametric==TRUE)	
 {
 
	cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "K     ", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   "); 
 
  }
 
 else 
 cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Ya    ", "Xa    " ,"Bo    ",  if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");  
  
 message(" ");	
 
 for(i in 1:length(half_life_values))
    {  
     for(k in 1:length(vy_values))
      { 
       if(half_life_values[i]==0) 
        {
         x.ols<-cbind(1, pred)
         beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
         vy <- vy_values[k];
        }
       else
        {
         a <- ln2/half_life_values[i];
         vy <- vy_values[k];
         x.ols<-cbind(1, pred)
         if(ultrametric==TRUE)
          beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
         else
         beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
        }

       ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
       con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
         if(half_life_values[i]==0)
          {
           a<-Inf	
           s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
           X<-cbind(1, pred)
           V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[2:(n.pred+1),])))                
           }
         else
         {
         if(ultrametric==TRUE)
          s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])) 
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
           X<-cbind(1, (1-(1-exp(-a*T))/(a*T))*pred)
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[2:(n.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))  

           V<-cm1+(s1*ta*cm2)+me.response+mv-mcov 
          }
         else
          {
           nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), (1-(1-exp(-a*T))/(a*T))*pred)
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[4:(n.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
           
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))   

           V<-cm1+(s1*ta*cm2)+me.response+mv-mcov
           
          }  
         } # END OF ELSE CONDITION FOR HALF-LIFE = 0

         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

         V.inverse<-solve(V)
         if(half_life_values[i]==0)
          {
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred+1))
           for(f in 1:(n.pred+1))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
         else
          {
         if(ultrametric==TRUE)
          {       
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred+1))
           for(f in 1:(n.pred+1))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
          } 
         else
          {
           beta.i<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)%*%(t(nu.X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred))
           for(f in 4:(n.pred+3))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
          }                          # END OF HALF-LIFE = 0 CONDITION #
        }                            # END OF ITERATED GLS REPEAT LOOP #
       beta1<-beta.i
   
      ### END OF ITERATED GLS ESTIMATION FOR BETA #
      
     if(half_life_values[i]==0)
      {
       s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))
       X<-cbind(1, pred)
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[2:(n.pred+1),])))
              V.inverse<-solve(V)
       eY<-X%*%beta1
       resid<-Y-eY;
       gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }
      else
      {
      if(ultrametric==TRUE)
       s1<-as.numeric(s.X%*%(beta1[2:(n.pred+1),]*beta1[2:(n.pred+1),]))   
      else
       s1<-as.numeric(s.X%*%(beta1[4:(n.pred+3),]*beta1[4:(n.pred+3),]))    
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
         }
       } 
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
         }
       }
      if(ultrametric==TRUE)
       {
        X<-cbind(1, (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[2:(n.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[2:(n.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))  

        V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
       }
      else
       {
        nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[4:(n.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))  
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[4:(n.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred))) 

        V<-cm1+(s1*ta*cm2)+me.response+mv-mcov 
       }
      V.inverse<-solve(V)
      if(ultrametric==TRUE)
       eY<-X%*%beta1
      else 
       eY<-nu.X%*%beta1
      resid<-Y-eY; 
      gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }  # END OF CONDITION FOR HALF-LIFE = 0 #
   
      print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4))) 
   
     }
    }
    
     
         
 x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
      }
    }
   gof=gof-ml
   
   

  # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #    
   
  if(alpha.est==Inf)
   {
     gls.beta1<-glsyx.beta1<- solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
     con.count<-0 # counter to break loop in the event of non-convergence
     repeat
      {
       s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))
       X<-cbind(1, pred)
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))
       V.inverse<-solve(V)
       beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(n.pred+1))
        for(f in 1:(n.pred+1))
         {
          if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
         } 
          if(sum(test)==0) break 
          con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
          gls.beta1<-glsyx.beta1<-beta.i
         }
      gls.beta1<-glsyx.beta1<-beta.i
      X<-cbind(1, pred)
      V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))
      pred.mean<-X%*%gls.beta1
      g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
      sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
      sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
      r.squared<-(sst-sse)/sst 
      
    }
    
  else
   {
   if(ultrametric==TRUE)
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
   else
    gls.beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y));
  con.count<-0;
    repeat
   {
    if(ultrametric==TRUE)
     s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])) 
    else
     s1<-as.numeric(s.X%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),]))   
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
     if(ultrametric==TRUE)
      {
       X<-cbind(1, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
       mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[2:(n.pred+1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))   
       mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[2:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred))) 

       V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
      } 
     else
      { 
       nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
       mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[4:(n.pred+3), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred))) 
        
       mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[4:length(gls.beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));   
      V<-cm1+(s1*ta*cm2)+me.response+mv-mcov 
      }
         
     V.inverse<-solve(V)
     
     
     if(ultrametric==TRUE)
      {       
       beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(n.pred+1))
       for(f in 1:(n.pred+1))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break  
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

          
         gls.beta1<-beta.i
       } 
      else
       {
        beta.i.var<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)
        beta.i<-beta.i.var%*%(t(nu.X)%*%V.inverse%*%Y)
        test<-matrix(nrow=(n.pred))
        for(f in 4:(n.pred+3))
         {
          if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
         } 
        if(sum(test)==0) break  
        con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
 
        beta1<-beta.i
       } 
   }       

        
         
                # END OF ITERATED GLS LOOP #


  # CODE FOR SST, SSE AND R-SQUARED #

  if(ultrametric==TRUE)
   gls.beta1<-beta.i 
  else
   {
    gls.beta1<-beta.i 
    ind.par<-matrix(data=0, nrow=N, ncol=4, dimnames=list(NULL, c("Bo", "Bi.Xia", "Yo", "Sum")))
    ind.par[,1]<-beta.i[1]*nu.X[,1]
    ind.par[,2]<-(beta.i[2]*nu.X[,2])
    ind.par[,3]<-beta.i[3]*nu.X[,3]
    ind.par[,4]<-ind.par[,1]+ind.par[,2]+ind.par[,3]
    mean.Bo=mean(ind.par[,4])
   }

  if(ultrametric==TRUE)
   {
    X<-cbind(1, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
    mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[2:(n.pred+1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))  
    mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[2:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))  
    V<-cm1+(s1*ta*cm2)+me.response+mv-mcov;
    pred.mean<-X%*%gls.beta1

   } 
   else
   { 
    nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred) 
    mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[4:(n.pred+3), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))  
    mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[4:length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))
    V<-cm1+(s1*ta*cm2)+me.response+mv-mcov 
    
    pred.mean<-nu.X%*%gls.beta1
   }

  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst  
  
  
  
  
                

  # FINAL EVOLUTIONARY REGRESSION USING BEST ALPHA AND VY ESTIMATES AND KNOWN VARIANCE MATRIX #
 
 
  if(ultrametric==TRUE)  s1<-as.numeric(s.X%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),]))
  else s1<-as.numeric(s.X%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),]));       
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
     
     if(ultrametric==TRUE)
    V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[2:(n.pred+1),]*gls.beta1[2:(n.pred+1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[2:length(gls.beta1),])))
    else
     V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[4:(n.pred+3),]*gls.beta1[4:(n.pred+3),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[4:length(gls.beta1),])));       
    
    X1<-cbind(1, pred) 
    V.inverse<-solve(V)
    ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
    ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)  
    glsyx.beta1<-ev.beta.i
  }                                         # END OF HALFLIFE 0 CONDITION #


  	
 } # END OF RANDOM COVARIATE REGRESSION ESTIMATION	
	
}# END OF FIXED COVARIATE, MIXED OR RANDOM MODELS PARAMETER ESTIMATION



 # EVALUATE IF IT IS A FIXED AND RANDOM COVARIATE ANCOVA OR REGRESSION MODEL ESTIMATE PARAMETERS WITH ITERATED GLS TO A) TAKE MEASUREMENT VARIANCE INTO ACCOUNT OR B) RANDOM EFFECTS INTO ACCOUNT IN THE CASE OF THE MIXED MODEL AND REGRESSION
    
if(model.type == "mmfANCOVA" || model.type=="mfReg") 
{
 # SET UP INITIAL MATRICES FOR MULTIPLE REGRESSION AND CALCULATE THETA AND SIGMA FOR RANDOM PREDICTOR / S
 
 pred<-data.frame(random.cov);
 n.pred<-length(pred[1,]);
 pred<-matrix(data=pred[!is.na(pred)], ncol=n.pred);
  if(is.null(me.random.cov)) me.pred<-matrix(data=0, nrow=N, ncol=n.pred) else me.pred<-matrix(data=me.random.cov[!is.na(me.random.cov)], ncol=n.pred);
 if(is.null(mecov.random.cov)) me.cov<-matrix(data=0, nrow=N, ncol=n.pred) else me.cov<-matrix(data=mecov.random.cov[!is.na(mecov.random.cov)], ncol=n.pred);
  
 s.X<-matrix(data=0, ncol=n.pred)  # PREDICTOR SIGMA
   for(i in 1:n.pred)
    {
    s.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[2]);
    }
      
   theta.X<-matrix(data=0, ncol=n.pred)  #PREDICTOR THETA
   for(i in 1:n.pred)
    {
     theta.X[,i] <- as.numeric(sigma.X.estimate(pred[,i],me.pred[,i], topology, times)[1]);
    }    
 
 # END OF RANDOM PREDICTOR THETA AND SIGMA ESTIMATES
 
 # FIXED COVARIATES
 
 fixed.pred<-data.frame(fixed.cov);
 n.fixed.pred<-length(fixed.pred[1,]);
 fixed.pred<-matrix(data=fixed.pred[!is.na(fixed.pred)], ncol=n.fixed.pred);
 if(is.null(me.fixed.cov)) me.fixed.pred<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.pred<- matrix(data=me.fixed.cov[!is.na(me.fixed.cov)], ncol=n.fixed.pred);
 if(is.null(mecov.fixed.cov)) me.fixed.cov<-matrix(data=0, nrow=N, ncol=n.fixed.pred) else me.fixed.cov<-matrix(data=me.cov.fixed.cov[!is.na(me.cov.fixed.cov)], ncol=n.fixed.pred);  
 
  
## INITIAL OLS ESTIMATES TO SEED ITERATED GLS

 if(model.type=="mfReg")  
  {
  	x.ols<-cbind(1, fixed.pred, pred);
    beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);  
    if(ultrametric == FALSE) beta1<-rbind(0, 0, beta1); # 2 additional parameter seeds for Ya and Xa    
  }
  
  if(model.type=="mmfANCOVA") 
  {
  	regime.specs<-fixed.fact;
  	n.fixed<-length(levels(as.factor(regime.specs)))
    regime.specs<-as.factor(regime.specs)
    x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred);  	beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);  
  }	
  
 
  
  # GRID ESTIMATION ROUTINE AND ITERATED GLS FOR MODELS THAT INCLUDE RANDOM EFFECTS
  
 if(model.type=="mmfANCOVA")   
 {
  
 cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", GS_head, if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   "); 
 
 
 message(" ");	
 	
 for(i in 1:length(half_life_values))
  {  
   for(k in 1:length(vy_values))
    { 
      if(half_life_values[i]==0) a<-1000000000000000000000 else a <- ln2/half_life_values[i];
      vy <- vy_values[k];
      X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred); 
       if(length(X[1,]) > length(beta1)) {beta1<-as.matrix(c(0, beta1)); n.fixed<-n.fixed+1}
       if(length(X[1,])< length(beta1)) {beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs))); print("The Ya parameter is dropped as its coefficient is too small");}        
                                  
    # CODE FOR ESTIMATING BETA USING ITERATED GLS 
    
       con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
        if(half_life_values[i]==0)
          {
           X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred);    
V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):(length(beta1)-n.pred),]*beta1[(n.fixed+1):(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),]))); 
}
        else
         {	
         	
          X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred);
          
                                              	
          s1<-as.numeric(s.X%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]));   
         
          
          for(p in 1:N)
           {
            for(q in 1:N)
             {
              if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
             }
           } 
          cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
          for(p in 1:N)
           {
            for(q in 1:N)
             {
             cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/ (a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
             }
           }

      mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1+n.fixed.pred):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));   
           mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));  
 mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(beta1[(n.fixed+1):(length(beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));  
           mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*beta1[(n.fixed+1):(length(beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));  

           
                      
           V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;     
          } # END OF If ELSE CONDITION FOR HALF-LIFE 0 OR NOT

         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #   
         
         V.inverse<-solve(V)
         if(half_life_values[i]==0)
          {
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
         else
          {
         
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(length(beta.i)))
           for(f in 1:(length(beta.i)))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
          }
         }  
                  
             
      ### END OF ITERATED GLS ESTIMATION FOR BETA #
      
     if(half_life_values[i]==0)
      {
       X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)   
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])))-diag(as.numeric(me.cov%*%(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[(n.fixed+1):(length(beta1)-n.pred),]*beta1[(n.fixed+1):(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[(n.fixed+1):(length(beta1)-n.pred),]))); 
       V.inverse<-solve(V)
       eY<-X%*%beta1
       resid<-Y-eY;
       gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }
      else
      {
       s1<-as.numeric(s.X%*%(beta1[(n.fixed+1+n.fixed.pred):length(beta1),]*beta1[(n.fixed+1+n.fixed.pred):length(beta1),])); 
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
         }
       } 
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
         }
       }
           X<-cbind(weight.matrix(a, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-a*T))/(a*T))*pred);          
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(n.fixed+1+n.fixed.pred):length(beta1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)));   
           mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(n.fixed+1+n.fixed.pred):length(beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)));  
 mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(beta1[(n.fixed+1):(length(beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));  
           mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*beta1[(n.fixed+1):(length(beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));  
           
           V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;
                   V.inverse<-solve(V)
     
       eY<-X%*%beta1
      
      resid<-Y-eY; 
      gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }  # END OF CONDITION FOR HALF-LIFE = 0 #
      print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4))) 
     }
    }
    
    
    
# END OF GRID SETUP,START OF GRID SEARCH FOR BEST ALPHA AND VY ESTIMATES #

   x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
      }
    }
   gof=gof-ml
   
  
   n.fixed<-length(levels(as.factor(regime.specs)))   ### reset before final regression
   
   
  # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #
   
  if(alpha.est==Inf || alpha.est >=1000000000000000000000)
   {
   	x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)   
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y) 
    con.count<-0;
    repeat
   	 {
   	 	
      s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]));
      X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.cov, intercept), pred)  
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]))); 
       V.inverse<-solve(V)
      beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(length(beta.i)))
       for(f in 1:(length(beta.i)))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break 
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
         gls.beta1<-beta.i
       }
     gls.beta1<-beta.i
     X<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)  
    V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),]))); 
     pred.mean<-X%*%gls.beta1
     g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
     sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
   
     sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
    
     r.squared<-(sst-sse)/sst 
         }
   else
   {
    x.ols<-cbind(weight.matrix(1000000000000000000000, topology, times, N, regime.specs, fixed.pred, intercept), pred)   
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y) 
    con.count<-0;
    
    X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred); 
       if(length(X[1,]) > length(gls.beta1)) {gls.beta1<-as.matrix(c(0, gls.beta1)); n.fixed<-n.fixed+1}
       if(length(X[1,])< length(gls.beta1)) {gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y);n.fixed<-length(levels(as.factor(regime.specs)))}        
   repeat
   {
   	
   	X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred);
     s1<-as.numeric(s.X%*%(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),]*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),])); 
   	     
 
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
    
                 
   mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)));   
           mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)));  
 mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));  
           mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));  
           
           V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;       
          V.inverse<-solve(V)
           beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(length(beta.i)))
       for(f in 1:(length(beta.i)))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break  
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

          
         gls.beta1<-beta.i
         
          X<-cbind(weight.matrix(alpha.est, topology, times, N, regime.specs, fixed.pred, intercept), (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
 
  mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)));   
           mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(n.fixed+1+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)));  
 mv.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.pred)*t(kronecker(gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred), ], rep(1, times=N))), ncol=n.fixed.pred)));  
           mcov.fixed<-diag(rowSums(matrix(data=as.numeric(me.fixed.cov)*t(kronecker(2*gls.beta1[(n.fixed+1):(length(gls.beta1)-n.pred),], rep(1, times=N))), ncol=n.fixed.pred)));  
           
           V<-cm1+(s1*ta*cm2)+me.response+mv+ mv.fixed-mcov-mcov.fixed;    
  pred.mean<-X%*%gls.beta1
  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst   
       }  
      }     

                # END OF ITERATED GLS LOOP #


            
                
      	
 } # END OF ESTIMATION MIXED MODEL ANCOVA	
 
 
 
if(model.type=="mfReg")   
 {
 
 if(ultrametric==TRUE)	
 {
 
	cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "K     ",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   "); 
 
  }
 
 else 
 cat(c("   ", "t1/2  ", "Vy    ", "Supp  ", "Ya    ", "Xa    " ,"Bo    ", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), sep="   ");  
  
 message(" ");	
 
 for(i in 1:length(half_life_values))
    {  
     for(k in 1:length(vy_values))
      { 
       if(half_life_values[i]==0) 
        {
         x.ols<-cbind(1, fixed.pred, pred)   
         beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
         vy <- vy_values[k];
        }
       else
        {
         a <- ln2/half_life_values[i];
         vy <- vy_values[k];
         x.ols<-cbind(1,fixed.pred, pred)
         if(ultrametric==TRUE)
          beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
         else
         beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y))
        }

       ### CODE FOR ESTIMATING BETA USING ITERATED GLS ###
       con.count<-0;  # Counter for loop break if Beta's dont converge #
       repeat
        {
         if(half_life_values[i]==0)
          {
           a<-Inf	
           s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
           X<-cbind(1, fixed.pred, pred)
           V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));   
             
           }
         else
         {
         if(ultrametric==TRUE)
          s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),])) 
         else
          s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),]*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),])) 
          
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
           X<-cbind(1, fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))

           V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));
          }
         else
          {
           nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
           mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
           
            mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+n.fixed.pred+3),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))  

           V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[4:(length(beta1)-n.pred),]*beta1[4:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));           
          }  
         } # END OF ELSE CONDITION FOR HALF-LIFE = 0

         # INTERMEDIATE ESTIMATION OF OPTIMAL REGRESSION #

         V.inverse<-solve(V)
         if(half_life_values[i]==0)
          {
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred+n.fixed.pred+1))
           for(f in 1:(n.pred+n.fixed.pred+1))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
         else
          {
         if(ultrametric==TRUE)
          {       
           beta.i<-pseudoinverse(t(X)%*%V.inverse%*%X)%*%(t(X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred+n.fixed.pred+1))
           for(f in 1:(n.pred+n.fixed.pred+1))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[f]=0 else test[f]=1
            } 
           if(sum(test)==0) break 
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
  
           beta1<-beta.i
          } 
         else
          {
           beta.i<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)%*%(t(nu.X)%*%V.inverse%*%Y)
           test<-matrix(nrow=(n.pred+n.fixed.pred))
           for(f in 4:(n.pred+n.fixed.pred+3))
            {
             if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
            } 
           if(sum(test)==0) break  
           con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

           beta1<-beta.i
          }
          }                          # END OF HALF-LIFE = 0 CONDITION #
        }                            # END OF ITERATED GLS REPEAT LOOP #
       beta1<-beta.i
   
      ### END OF ITERATED GLS ESTIMATION FOR BETA #
      
     if(half_life_values[i]==0)
      {
       s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))
       X<-cbind(1, fixed.pred,pred)
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred++1),]*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),])))-diag(as.numeric(me.cov%*%(2*beta1[(2+n.fixed.pred):(n.pred+n.fixed.pred+1),]))) + diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),]))) 
              V.inverse<-solve(V)
       eY<-X%*%beta1
       resid<-Y-eY;
       gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }
      else
      {
      if(ultrametric==TRUE)
       s1<-as.numeric(s.X%*%(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))   
      else
       s1<-as.numeric(s.X%*%(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]))    
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-a*ta[q,p]))/(a*ta[q,p]);
         }
       } 
      cm1<-(s1/(2*a)+vy)*(1-exp(-2*a*ta))*exp(-a*tij);
      for(p in 1:N)
       {
        for(q in 1:N)
         {
          cm2[p,q]<-(((1-exp(-a*T[p]))/(a*T[p]))*((1-exp(-a*T[q]))/(a*T[q]))-(exp(-a*tia[p, q])*(1-exp(-a*T[p]))/(a*T[q])+ exp(-a*tja[p, q])*(1-exp(-a*T[p]))/(a*T[p]))*(num.prob[p,q]));
         }
       }
      if(ultrametric==TRUE)
       {
        X<-cbind(1, fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred))) 
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred))) 
        V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(beta1)-n.pred),]*beta1[2:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[2:(length(beta1)-n.pred),])));
        }
      else
       {
        nu.X<-cbind(1-exp(-a*T), 1-exp(-a*T)-(1-(1-exp(-a*T))/(a*T)), exp(-a*T), fixed.pred, (1-(1-exp(-a*T))/(a*T))*pred)
        mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-a*T))/(a*T)))^2), ncol=n.pred)))  
        mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred))) 

        V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(beta1[4:(length(beta1)-n.pred),]*beta1[4:(length(beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*beta1[4:(length(beta1)-n.pred),])));       
        }  
      V.inverse<-solve(V)
      if(ultrametric==TRUE)
       eY<-X%*%beta1
      else 
       eY<-nu.X%*%beta1
      resid<-Y-eY; 
      gof[i, k] <- -N/2*log(2*pi)-0.5*log(det(V))-0.5*(t(resid) %*% V.inverse%*%resid); 
      }  # END OF CONDITION FOR HALF-LIFE = 0 #
   
      print(as.numeric(round(cbind(if(a!=0)log(2)/a else 0.00, vy, gof[i,k], t(beta1)), 4))) 
   
     }
    }
    
     
         
 x<-rev(half_life_values)
   y<-vy_values   
   z<-gof;
   ml<-max(z);
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]==ml){alpha.est=log(2)/half_life_values[i]; vy.est=vy_values[j]}
      }
    }
   for(i in 1:length(half_life_values))
    {
     for(j in 1:length(vy_values))
      {
       if(gof[i,j]<=ml-support)gof[i, j]=ml-support;
      }
    }
   gof=gof-ml
   
   

  # FINAL OPTIMAL REGRESSION USING BEST ALPHA AND VY ESTIMATES #    
   
  if(alpha.est==Inf)
   {
     gls.beta1<-glsyx.beta1<- solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
     con.count<-0 # counter to break loop in the event of non-convergence
     repeat
      {
       s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
       X<-cbind(1, fixed.pred, pred)
       V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),]))); 
       V.inverse<-solve(V)
       beta.i.var<-ev.beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(n.pred+n.fixed.pred+1))
        for(f in 1:(n.pred+1+n.fixed.pred))
         {
          if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
         } 
          if(sum(test)==0) break 
          con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
          gls.beta1<-glsyx.beta1<-beta.i
         }
      gls.beta1<-glsyx.beta1<-beta.i
      X<-cbind(1, fixed.pred,pred)
      V<-diag(rep(vy, times=N))+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),]))); 
      pred.mean<-X%*%gls.beta1
      g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
      sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
      sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
      r.squared<-(sst-sse)/sst 
      
    }
    
  else
   {
   if(ultrametric==TRUE)
    gls.beta1<-solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y)
   else
    gls.beta1<-rbind(0, 0, solve(t(x.ols)%*%x.ols)%*%(t(x.ols)%*%Y));
  con.count<-0;
    repeat
   {
    if(ultrametric==TRUE)
     s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])) 
    else
     s1<-as.numeric(s.X%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]))   
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
     if(ultrametric==TRUE)
      {
       X<-cbind(1, fixed.pred,(1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
       mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))   
       mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred))) 

       V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])));      } 
     else
      { 
       nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T), fixed.pred, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
       mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred))) 
        
       mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(4+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-a*T))/(a*T)))), ncol=n.pred)))
      V<-cm1+(s1*ta*cm2)+me.response+mv-mcov + diag(as.numeric(me.fixed.pred%*%(gls.beta1[4:(length(gls.beta1)-n.pred),]*gls.beta1[4:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])));
       }
         
     V.inverse<-solve(V)
     
     
     if(ultrametric==TRUE)
      {       
       beta.i.var<-pseudoinverse(t(X)%*%V.inverse%*%X)
       beta.i<-beta.i.var%*%(t(X)%*%V.inverse%*%Y)
       test<-matrix(nrow=(n.pred+1+n.fixed.pred))
       for(f in 1:(n.pred+1+n.fixed.pred))
        {
         if(abs(as.numeric(beta.i[f]-gls.beta1[f]))<=convergence) test[f]=0 else test[f]=1
        } 
         if(sum(test)==0) break  
         con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  

          
         gls.beta1<-beta.i
       } 
      else
       {
        beta.i.var<-pseudoinverse(t(nu.X)%*%V.inverse%*%nu.X)
        beta.i<-beta.i.var%*%(t(nu.X)%*%V.inverse%*%Y)
        test<-matrix(nrow=(n.pred))
        for(f in 4:(n.pred+3+n.fixed.pred))  
         {
          if(abs(as.numeric(beta.i[f]-beta1[f]))<=convergence) test[(f-3)]=0 else test[(f-3)]=1
         } 
        if(sum(test)==0) break  
        con.count=con.count+1
           if(con.count >= 50) 
            {
             message("Warning, Beta estimates did not converge after 50 iterations, last estimates printed out")
             break
            }  
 
        beta1<-beta.i
       } 
   }       

        
         
                # END OF ITERATED GLS LOOP #


  # CODE FOR SST, SSE AND R-SQUARED #

  if(ultrametric==TRUE)
   gls.beta1<-beta.i 
  else
   {
    gls.beta1<-beta.i 
    ind.par<-matrix(data=0, nrow=N, ncol=4, dimnames=list(NULL, c("Bo", "Bi.Xia", "Yo", "Sum")))
    ind.par[,1]<-beta.i[1]*nu.X[,1]
    ind.par[,2]<-(beta.i[2]*nu.X[,2])
    ind.par[,3]<-beta.i[3]*nu.X[,3]
    ind.par[,4]<-ind.par[,1]+ind.par[,2]+ind.par[,3]
    mean.Bo=mean(ind.par[,4])
   }

  if(ultrametric==TRUE)
   {
    X<-cbind(1, fixed.pred,(1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred)
    mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))  
    mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred)))  
    V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))
    pred.mean<-X%*%gls.beta1

   } 
   else
   { 
    nu.X<-cbind(1-exp(-alpha.est*T), 1-exp(-alpha.est*T)-(1-(1-exp(-alpha.est*T))/(alpha.est*T)), exp(-alpha.est*T),fixed.pred, (1-(1-exp(-alpha.est*T))/(alpha.est*T))*pred) 
    mv<-diag(rowSums(matrix(data=as.numeric(me.pred)*t(kronecker(beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred), ], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))^2), ncol=n.pred)))  
    mcov<-diag(rowSums(matrix(data=as.numeric(me.cov)*t(kronecker(2*beta1[(4+n.fixed.pred):length(beta1),], (1-(1-exp(-alpha.est*T))/(alpha.est*T)))), ncol=n.pred))) 
    V<-cm1+(s1*ta*cm2)+me.response+mv-mcov+ diag(as.numeric(me.fixed.pred%*%(gls.beta1[4:(length(gls.beta1)-n.pred),]*gls.beta1[4:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])))
    
    pred.mean<-nu.X%*%gls.beta1
   }

  g.mean<-(t(rep(1, times=N))%*%solve(V)%*%Y)/sum(solve(V));
  sst<-t(Y-g.mean)%*% solve(V)%*%(Y-g.mean)
  sse<-t(Y-pred.mean)%*%solve(V)%*%(Y-pred.mean)
  r.squared<-(sst-sse)/sst  
  
  
  
  
                

  # FINAL EVOLUTIONARY REGRESSION USING BEST ALPHA AND VY ESTIMATES AND KNOWN VARIANCE MATRIX #
 
 
  if(ultrametric==TRUE)  s1<-as.numeric(s.X%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]))
  else s1<-as.numeric(s.X%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]));       
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        if(ta[q,p]==0)num.prob[q,p]=1 else num.prob[q,p]=(1-exp(-alpha.est*ta[q,p]))/(alpha.est*ta[q,p])
       }
     } 
    cm1<-(s1/(2*alpha.est)+vy.est)*(1-exp(-2*alpha.est*ta))*exp(-alpha.est*tij)
    for(p in 1:N)
     {
      for(q in 1:N)
       {
        cm2[p,q]<-(((1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*((1-exp(-alpha.est*T[q]))/(alpha.est*T[q]))-(exp(-alpha.est*tia[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[q])+ exp(-alpha.est*tja[p, q])*(1-exp(-alpha.est*T[p]))/(alpha.est*T[p]))*(num.prob[p,q]))
       }
     }
     
     if(ultrametric==TRUE)
    V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),]*gls.beta1[(2+n.fixed.pred):(n.pred+1+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(2+n.fixed.pred):length(gls.beta1),]))) + diag(as.numeric(me.fixed.pred%*%(gls.beta1[2:(length(gls.beta1)-n.pred),]*gls.beta1[2:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[2:(length(gls.beta1)-n.pred),])))
        else
     V<-cm1+(s1*ta*cm2)+me.response+diag(as.numeric(me.pred%*%(gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),]*gls.beta1[(4+n.fixed.pred):(n.pred+3+n.fixed.pred),])))-diag(as.numeric(me.cov%*%(2*gls.beta1[(4+n.fixed.pred):length(gls.beta1),])))+ diag(as.numeric(me.fixed.pred%*%(gls.beta1[4:(length(gls.beta1)-n.pred),]*gls.beta1[4:(length(gls.beta1)-n.pred),])))-diag(as.numeric(me.fixed.cov%*%(2*gls.beta1[4:(length(gls.beta1)-n.pred),])))
    
    X1<-cbind(1, fixed.pred, pred) 
    V.inverse<-solve(V)
    ev.beta.i.var<-pseudoinverse(t(X1)%*%V.inverse%*%X1)
    ev.beta.i<-ev.beta.i.var%*%(t(X1)%*%V.inverse%*%Y)  
    glsyx.beta1<-ev.beta.i
  }                                         # END OF HALFLIFE 0 CONDITION #


  	
 } # END OF RANDOM AND FIXED COVARIATE REGRESSION ESTIMATION	
	
}# END OF FIXED AND RANDOM COVARIATE ANCOVA AND REGRESSION PARAMETER ESTIMATION
#### END OF NEW CODE


# PLOT THE SUPPORT SURFACE FOR HALF-LIVES AND VY


   if(length(half_life_values) > 1 && length(vy_values) > 1){
    z1<-gof
    for(i in 1:length(vy_values)){
     h.lives[,i]=rev(z1[,i])
     }
     z<-h.lives
    op <- par(bg = "white")
    persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "NA")
     persp(x, y, z, theta = 30, phi = 30, expand = 0.5, col = "NA",
            ltheta = 120, shade = 0.75, ticktype = "detailed",
            xlab = "half-life", ylab = "vy", zlab = "log-likelihood") -> res
            }
            
# MODEL OUTPUT

# alpha, half-lives, correction factor, v


  message("==================================================")
   half.life<-log(2)/alpha.est
   c.factor<-mean(1-(1-exp(-alpha.est*T))/(alpha.est*T))
   modeloutput<-matrix(data=0, nrow=4, ncol=1, dimnames=list(c("Rate of adaptation ", "Phylogenetic half-life ","Phylogenetic correction factor", "Stationary variance "), "    Estimate"))
   modeloutput[1, 1]=alpha.est; modeloutput[2, 1]=half.life; modeloutput[3,1]=c.factor; modeloutput[4,1]=vy.est;   ##### Rememeber to output s.X 





 modfit<-matrix(data=0, nrow=7, ncol=1, dimnames=list(c("Support", "AIC", "AICc", "SIC", "r squared", "SST", "SSE"),("Value")))
 
 
   #if(ultrametric==TRUE) n.par=1+n.pred else n.par=3+n.pred
   
  if(model.type=="ffANOVA" || model.type=="fReg" || model.type=="ffANCOVA") n.par<-length(gls.beta0)
  if(model.type == "mmANCOVA" || model.type=="rReg" || model.type=="mfReg" || model.type=="mmfANCOVA")   n.par<-length(beta1)
   
   modfit[1,1]=ml 
   modfit[2,1]=-2*ml+2*(2+n.par)
   modfit[3,1]=modfit[2,1]+(2*(2+n.par)*((2+n.par)+1))/(N-(2+n.par)-1)
   modfit[4,1]=-2*ml+log(N)*(2+n.par)
   modfit[5,1]=r.squared*100 
   modfit[6,1]=sst
   modfit[7,1]=sse
   
   message("");
   message("BEST ESTIMATES & MODEL FIT");message("");  
   message("==================================================");
   message("MODEL PARAMETERS"); 
   print(modeloutput);message("");
   
# predictor means and variances for random predictors   

if(model.type == "mmANCOVA" || model.type=="rReg" || model.type=="mfReg" || model.type=="mmfANCOVA")
{
print(matrix(data=rbind(theta.X, s.X), nrow=2, ncol=n.pred, dimnames=list(c("Predictor theta", "Predictor variance"), if(n.pred==1) deparse(substitute(random.cov)) else colnames(random.cov))));
message("");
}

# PRIMARY OPTIMA OR REGRESSION SLOPE ESTIMATES

 message("--------------------------------------------------");
 message("PRIMARY OPTIMA");message(""); 

 
if(model.type=="IntcptReg")
{
if(ultrametric==TRUE || alpha.est==Inf || alpha.est>=1000000000000000){
   Intercept<-matrix(nrow=1, ncol=2, dimnames=list(("Theta_global"), c("Estimate", "Std.error")))
   Intercept[,1]<-gls.beta0
   Intercept[,2]<-sqrt(beta.i.var)}
   else {
   Intercept<-matrix(data=0, nrow=2, ncol=1, dimnames=list(c("Bo", "Ya"), ("     Estimate")))
   Intercept[1,1]<-beta.i[1]
   Intercept[2,1]<-beta.i[2]
   }
    print(Intercept); message("")  
 }    
   
if(model.type=="ffANOVA")
{
 std<-sqrt(diag(beta.i.var))

 optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(colnames(X), c("Estimates", "Std.error")));
 optima[,1] = gls.beta0;
 optima[,2] = std;
 
  
 reg <- set.of.regimes(topology,regime.specs);
 root.reg<-as.character(regime.specs[times==0])
 nonroot.reg<-as.character(reg[reg != root.reg])
   
    
 if(is.null(intercept)) 
    {
     if(ncol(X) == length(reg)) message ("The ancestral state (Ya) parameter was dropped from this model as there is not enough information to estimate it")  else
     if(ncol(X)<length(reg)) message ("Ya and the parameter at the root were dropped") else  
    message("this model does not drop Ya as it may influence the other parameters")
    }
    else
     {
      if(intercept=="root") message(root.reg, " ", "mapped to the root of the tree and includes the coefficent for the ancestral state (Ya)") else
      message("you set the intercept coefficent to a value of", " ", intercept,". Ya is not the true ancestral state anymore")
     }
     print(optima);message("");
 }  
 
 
if(model.type== "fReg")
{
 std<-sqrt(diag(beta.i.var))
 
 optima<-matrix(data=0, nrow=(nrow(gls.beta0)), ncol=2, dimnames=list(c("Bo", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Estimate", "Std. Error")))
 optima[,1] = gls.beta0;
 optima[,2] = std;
    
 print(optima);message("");
}  

if(model.type=="ffANCOVA")
{
 std<-sqrt(diag(beta.i.var))
 
 
 optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov)), c("Estimates", "Std.error")));
  optima[,1] = gls.beta0;
 optima[,2] = std;
    
 print(optima);message("");
}  

 
 
if(model.type  == "mmANCOVA")
{
std<-sqrt(diag(beta.i.var))	

if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")))
else
	
optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")));

 optima[,1] = gls.beta1;
 optima[,2] = std;
 print(optima)
 	
}

if(model.type  == "mmfANCOVA")
{
std<-sqrt(diag(beta.i.var))	

if(length(X[1,]) > length(x.ols[1,])) optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(c("Ya",as.character(levels(fixed.fact))),if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")))
else
	
optima<-matrix(data=0, nrow=ncol(X), ncol=2, dimnames = list(c(as.character(levels(fixed.fact)), if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimates", "Std.error")));

 optima[,1] = gls.beta1;
 optima[,2] = std;
 print(optima)
 	
}

 

if(model.type=="rReg")   
{
if(ultrametric==TRUE || alpha.est == Inf)	
 opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
 else
 {
  if(alpha.est != Inf)  
 opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error"))) 
 else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))}
 
     opreg[,1] =round(gls.beta1, 5)
     opreg[,2]= round(sqrt(diag(beta.i.var)),5)
     
 if(model.type=="rReg")
 {
 	
 	evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept", if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
 
     
     evreg[,1] =round(glsyx.beta1, 5)
     evreg[,2]= round(sqrt(diag(ev.beta.i.var)),5)    
      
     message("Evolutionary regression"); message("")   
     print(evreg);
     message("");
 }    
     message("Optimal regression"); message("")
     print(opreg);
     
if(model.type=="rReg" && ultrametric==TRUE && alpha.est != Inf)
{
     message("")
     message("Decomposition of K assuming Ya = Xa to get the optimal regression intercept Bo")
     message("")

  bo<-opreg[1,1] + (c.factor-1)*(sum(gls.beta1[-1]*theta.X)) 
  print(bo)
  message("")
  message("(Use this as the intercept when plotting the regression line)")
  
  message("")
}  
} 


if(model.type=="mfReg")   
{
if(ultrametric==TRUE || alpha.est == Inf)	
 opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
 else
 {
  if(alpha.est != Inf)  
 opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("Xa", "Bo","Ya" ,if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error"))) 
 else opreg<-matrix(data=0, nrow=(nrow(gls.beta1)), ncol=2, dimnames=list(c("K", if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov),if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))}
 
     opreg[,1] =round(gls.beta1, 5)
     opreg[,2]= round(sqrt(diag(beta.i.var)),5)
     
 if(model.type=="mfReg")
 {
 	
 	evreg<-matrix(data=0, nrow=(nrow(glsyx.beta1)), ncol=2, dimnames=list(c("Intercept",if(is.null(dim(fixed.cov))) deparse(substitute(fixed.cov)) else colnames(fixed.cov), if(is.null(dim(random.cov))) deparse(substitute(random.cov)) else colnames(random.cov)), c("Estimate", "Std. Error")))
 
     
     evreg[,1] =round(glsyx.beta1, 5)
     evreg[,2]= round(sqrt(diag(ev.beta.i.var)),5)    
      
     message("Evolutionary regression"); message("")   
     print(evreg);
     message("");
 }    
     message("Optimal regression"); message("")
     print(opreg);
     
if(model.type=="mfReg" && ultrametric==TRUE && alpha.est != Inf)
{
     message("")
     message("Decomposition of K assuming Ya = Xa to get the optimal regression intercept Bo")
     message("")

  bo<-opreg[1,1] + (c.factor-1)*(sum(gls.beta1[-(1:(1+n.fixed.pred))]*theta.X)) 
  print(bo)
  message("")
  message("(Use this as the intercept when plotting the regression line)")
  
  message("")
}  
} 
 
 message("--------------------------------------------------"); 
   message("MODEL FIT");message("");
   print(modfit); message("");
   message("=================================================="); 

 
} # END OF MODEL FITTING FUNCTION

