`parse.tree` <-
function(topology, times){
 term<-terminal.twigs(topology);
 N <- length(term);
 dm<-distance.matrix(topology, times);
 bt<-branch.times(topology, times);
 pt<-list(N=N, term=term, dm=dm, bt=bt);
 return(pt);
}

