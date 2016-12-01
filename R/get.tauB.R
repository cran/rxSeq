`get.tauB` = 
function(i, n, n0B, geneids, min.cnt=50, exclude.prop=.05, filt="ENSMUSG00000086503"){
    if(!is.null(filt)){
      filter = which(!is.na(match(geneids, filt)))
    }else{
      filter = integer(0)
    }
    if(length(filter) >= nrow(n)){
      stop("all genes were excluded by the Xist.ID filter")
    }    
    if(length(filter)>0){
      n.f = n[ -filter, i]
      sub.n0B = n0B[ -filter, i]
    }else{
      n.f = n[, i]
      sub.n0B = n0B[, i]
    }
    sub.n0B = sub.n0B[n.f >= min.cnt]
    n.f = n.f[n.f >= min.cnt]
    sub.n0A = n.f - sub.n0B
    flag1 = sub.n0A/n.f >= exclude.prop
    flag2 = sub.n0B/n.f >= exclude.prop    
    flag = flag1  &  flag2
    if(sum(flag) == 0){
      stop("too stringent criterias for the given counts")
    }
    counts = sum(sub.n0B[flag])
    counts2 = sum(n.f[flag])
    out = c(median(sub.n0B[flag]/n.f[flag]), counts/counts2, length(sub.n0B), sum(flag))

    return(out)
}
