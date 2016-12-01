`get.tausB` = 
function(n, n0B, geneids, min.cnt=50, exclude.prop=.05, Xist.ID="ENSMUSG00000086503"){
  out = matrix(sapply(1:ncol(n), get.tauB, n=n, n0B=n0B, geneids=geneids, min.cnt=min.cnt, exclude.prop=exclude.prop, filt=Xist.ID), nrow=4)
  if(any(unlist(sapply(out, is.null)))){
    stop("too stringent criterias for the given counts")
  }
  colnames(out) = colnames(n0B)
  rownames(out) = c("med.tauB", "ave.tauB", "all.genes", "used.genes")
  return(out)
}
