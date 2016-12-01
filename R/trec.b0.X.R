`trec.b0.X` = 
function(yi, ind.lst, X, twosex, l.tau.r, l.tau.a, 
               start, iphi=1, theta=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  if(twosex){zeroes = c(0, 0)}else{zeroes = 0}
  par0 = start
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      iphi_i = iphi
      log.lik0 = ll.toRCI.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
      for(i in 1:maxiter){
        out = optim(par0, ll.toRCI.b0.X, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
        par0 = out$par
        log.lik1 = ll.toRCI.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
        older = c(zeroes, par0, log.lik1, iphi_i)        
        iphi_i = optimize(f=ll.toRCI.iphi.X, c(tech.ctrl$iphi_l, tech.ctrl$iphi_u), 
              yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, betas=c(zeroes, par0), l.tau.r=l.tau.r, l.tau.a=l.tau.a)$minimum    
        log.lik1 = ll.toRCI.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
        newer = c(zeroes, par0, log.lik1, iphi_i)      
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
