`nLogLik.trecase.A` = function(coef, phi, theta, rc, genei, hessian=FALSE){
#coef, phi, theta
  index = rc$index
  y = rc$y[genei,]
  n = rc$n[genei,]
  n0B = rc$n0B[genei,]
  kappas = rc$kappas

  out = input.checks.A(coef=coef, phi=phi, index=index, y=y, kappas=kappas, theta=theta, n=n, n0B=n0B, trecase=TRUE)
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
  xs = index[1:length(n)];xs[xs == 5] = 1;xs[xs %in% c(2, 6)] = -1
  if(twosex){
    sex.tot = intercept;sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dom*sex.tot)
    sex.ase = sex.tot[1:length(n)]
  }else{
    X = cbind(intercept, kappas, dom)
    sex.ase = NULL
  }
  nll = ll.jRCI.A(coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase,  
            ni=n, ni0=n0B, xs=xs, iphi=1/phi, theta=theta)

  if(hessian){
    tag = tryCatch({    
      hess = hessian(ll.jRCI.A, x=coef, yi=y, ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase,  
            ni=n, ni0=n0B, xs=xs, iphi=1/phi, theta=theta)
      0
    }, error=function(e) {
      warning("estimation of hessian failed")
      1
    })    
  }
  return(list(nll=nll, hess=hess))
}


