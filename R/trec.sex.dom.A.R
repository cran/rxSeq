`trec.sex.dom.A` = 
function(yi, ind.lst, X, start, iphi=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  par0 = start
  seq = 1:ncol(X)
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      res = trec.full.A(yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, start=par0, iphi=iphi, maxiter=maxiter, eps=eps, tech.ctrl=tech.ctrl)
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