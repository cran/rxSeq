`nLogLik.trec.X` = function(coef, phi, rc, genei, hessian=FALSE){
  index = rc$index
  y = rc$y[genei,]
  kappas = rc$kappas
  tausB = rc$tausB
  gene.switch = rc$geneids[genei] %in% rc$genes.switch  

  out = input.checks.X(coef=coef, phi=phi, index=index, y=y, kappas=kappas, trecase=FALSE)
  hess = NULL
  twosex = max(index)>4
  if(twosex){
    ind.lst = lapply(1:8, function(i){which(index %in% i)})
  }else{
    ind.lst = lapply(1:4, function(i){which(index %in% i)})    
  }
  if(gene.switch){
    l.tau.a = log(tausB)
    l.tau.b = log(1 - tausB)
  }else{
    l.tau.a = log(1 - tausB)
    l.tau.b = log(tausB)
  }
  l.tau.r = l.tau.b - l.tau.a  

  nsamples = length(index)
  intercept = poo = dev = dom = rep(0, length(index))
  intercept = intercept + 1
  dom[index %in% c(1, 2)] = 1
  dev[index %in% c(5, 6)] = 1  
  poo[index == 1] = 1
  poo[index == 2] = -1
  if(twosex){
    sex.tot = intercept; sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dev, poo)
  }else{
    X = cbind(intercept, kappas, dom, poo)
  }
  nll = ll.toRCI.X(coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, 
         iphi=1/phi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)  

  if(hessian){
    tag = tryCatch({    
      hess = hessian(ll.toRCI.X, x=coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, 
         iphi=1/phi, l.tau.r=l.tau.r, l.tau.a=l.tau.a)  
      0
    }, error=function(e) {
      warning("estimation of hessian failed")
      1
    })    
  }
  return(list(nll=nll, hess=hess))
}


