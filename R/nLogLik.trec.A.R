`nLogLik.trec.A` = function(coef, phi, rc, genei, hessian=FALSE){
  index = rc$index
  y = rc$y[genei,]
  kappas = rc$kappas

  out = input.checks.A(coef=coef, phi=phi, index=index, y=y, kappas=kappas, trecase=FALSE)
  hess = NULL
  twosex = max(index)>4
  if(twosex){
    ind.lst = lapply(1:8, function(i){which(index %in% i)})
  }else{
    ind.lst = lapply(1:4, function(i){which(index %in% i)})    
  }
  nsamples = length(index)
  intercept = rep(1, nsamples)
  dom = as.numeric(index %in% c(1, 2, 5, 6))
  poo = rep(0, length(index))
  poo[index %in% c(1, 5)] = 1
  poo[index %in% c(2, 6)] = -1
  if(twosex){
    sex.tot = intercept;sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dom*sex.tot, poo, poo*sex.tot)
  }else{
    X = cbind(intercept, kappas, dom, poo)
  }
  nll = ll.toRCI.A(coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, iphi=1/phi)

  if(hessian){
    tag = tryCatch({    
      hess = hessian(ll.toRCI.A, x=coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, iphi=1/phi)
      0
    }, error=function(e) {
      warning("estimation of hessian failed")
      1
    })    
  }
  return(list(nll=nll, hess=hess))
}


