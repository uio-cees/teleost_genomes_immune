`tree.layout` <-
function (topology) {
  root <- which(topology==0);
  return(arrange.tree(root,topology));
}

