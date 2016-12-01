`trec.b0.A` = 
function(yi, ind.lst, X, twosex, 
               start, iphi=1, theta=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  if(twosex){zeroes = c(0, 0)}else{zeroes = 0}
  par0 = start
  iphi_l = tech.ctrl$iphi_l
  repeat{
    triali = triali + 1
    tag = tryCatch({    
      iphi_i = iphi
      log.lik0 = ll.toRCI.b0.A(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i)
      for(i in 1:maxiter){
        out = optim(par0, ll.toRCI.b0.A, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i)
        par0 = out$par
        log.lik1 = ll.toRCI.b0.A(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i)
        older = c(zeroes, par0, log.lik1, iphi_i)        
        iphi_i = optimize(f=ll.toRCI.iphi.A, c(tech.ctrl$iphi_l, tech.ctrl$iphi_u), 
                       yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, betas=c(zeroes, par0))$minimum    
        log.lik1 = ll.toRCI.b0.A(par0, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi_i)
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
#    if((tag == 0) | (triali >= tech.ctrl$maxtrial))break;
    if((tag == 0) | (triali >= tech.ctrl$maxtrial))break;
    par0 = rnorm(length(par0), start, 1)
  }
  if(tag == 1){
    ret = NULL
  }
  return(ret)
}
