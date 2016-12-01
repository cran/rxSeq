`trecase.dom.X` = 
function(yi, ind.lst, X, twosex, ni, ni0, xs, l.tau.r, l.tau.a, 
           start, iphi=1, theta=1, maxiter=100, eps=1E-3, tech.ctrl){
  ret = NULL
  triali = 0
  par0 = start
  if(twosex){
    seq = 1:(3 + ncol(X))
  }else{
    seq = 1:(2 + ncol(X))
  }
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      res = trecase.full.X(yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, ni=ni, ni0=ni0, xs=xs, l.tau.r=l.tau.r, l.tau.a=l.tau.a,  
                     start=par0, iphi=iphi, theta=theta, maxiter=maxiter, eps=eps, tech.ctrl=tech.ctrl)
      0
    }, error=function(e) {
      1
    })
    if((tag == 0) | (triali >= tech.ctrl$maxtrial))break;
    par0 = rnorm(length(start), start, 1)
  }
  if(tag == 1){
    ret = NULL
  }else{
    if(!is.null(res)){
      ret = c(res[seq], 0, res[-seq])
    }
  }
  return(ret)      
}
