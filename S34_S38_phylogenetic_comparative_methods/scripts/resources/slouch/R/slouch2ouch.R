`slouch2ouch` <-
function(nodes, ancestors, times, labels)
 {
  ancestors[1]<-NA		
  tree<-ouchtree(nodes, ancestors, times, labels)	
  return(tree)
 }

