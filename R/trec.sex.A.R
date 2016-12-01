`trec.sex.A` = 
function(yi, ind.lst, X, start, iphi=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  par0 = start
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      iphi_i = iphi
      log.lik0 = ll.toRCI.sex.A(par0=par0, yi=yi, ind.lst=ind.lst, X=X, iphi=iphi_i)
      for(i in 1:maxiter){
        out = optim(par0, ll.toRCI.sex.A, yi=yi, ind.lst=ind.lst, X=X, iphi=iphi_i)
        log.lik1 = ll.toRCI.sex.A(par0=par0, yi=yi, ind.lst=ind.lst, X=X, iphi=iphi_i)    
        older = c(par0[1], par0[1:3], 0, par0[4], 0, par0[5], 0, log.lik1, iphi_i)
        par0 = out$par
        iphi_i = optimize(f=ll.toRCI.iphi.A, c(tech.ctrl$iphi_l, tech.ctrl$iphi_u), 
                       yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, betas=c(par0[1], par0[1:5]))$minimum
        log.lik1 = ll.toRCI.sex.A(par0, yi=yi, ind.lst=ind.lst, X=X, iphi=iphi_i)    
        newer = c(par0[1], par0[1:3], 0, par0[4], 0, par0[5], 0, log.lik1, iphi_i)      
        if((log.lik0 - log.lik1)<eps)break;
        log.lik0 = log.lik1
      }
      if(log.lik0 < log.lik1){
        ret = older
      }else{
        ret = newer
      }
      0
    }, error=function(e) {
      1
    })
    if((tag == 0) | (triali >= tech.ctrl$maxtrial))break;
    par0 = rnorm(length(par0), start, 1)
  }
  if(tag == 1){
    ret = NULL
  }
  return(ret)
}
