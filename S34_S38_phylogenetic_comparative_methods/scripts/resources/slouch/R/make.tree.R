`make.tree` <-
function(n.tips, shape=NULL, stretch=NULL, ouch=FALSE){
 nodes<-n.tips*2-1
 int.nodes<-nodes-n.tips
 species<-c(rep(NA, times=int.nodes), colours()[1:n.tips])
 if(is.null(shape))
 anc<-sort(c(0, seq(1:int.nodes), seq(1:int.nodes)))
 else
 anc<-c(0,seq(1:int.nodes), rev(seq(1:int.nodes)))
 times<-c(0, anc[2:int.nodes], rep(anc[int.nodes+1], times=n.tips))
 times<-times/max(times)
 if(!is.null(stretch)) times<-times^stretch
 
 # decide to output ouch or slouch tree
 
 if(ouch==TRUE)
 {
 anc[1]<-NA	
 tree<-ouchtree(nodes=1:nodes, ancestors=anc, times=times,labels=species)
 }
 else
 {
 tree<-data.frame(node=1:nodes, ancestor=anc, time=times, species=species)
 }
 return(tree)
 }

