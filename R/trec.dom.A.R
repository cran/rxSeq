`trec.dom.A` = 
function(yi, ind.lst, X, twosex, start, iphi=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  seq = 1:ncol(X)
  par0 = start
  if(twosex){
    zeroes = c(0, 0)
  }else{
    zeroes = 0
  }
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      res = trec.full.A(yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, start=par0, 
                        iphi=iphi, maxiter=maxiter, eps=eps, tech.ctrl=tech.ctrl)
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
      ret = c(res[seq], zeroes, res[-seq]) 
    }
  }
  return(ret)  
}
