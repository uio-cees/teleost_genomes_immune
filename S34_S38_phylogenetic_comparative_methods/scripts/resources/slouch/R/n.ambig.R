`n.ambig` <-
function(pars.states)
{
 states<-pars.states$Final.states
 N<-length(states)
  count<-0
  ambig<-NA
   for(i in 1:N)
   {
  if(length(states[[i]]) >=2) 
  {
  count=count+1;
  ambig<-c(ambig, i)
  }
  
   }
    return(list(N.ambiguous=count, ambiguous.nodes=ambig[-1]))	
}

