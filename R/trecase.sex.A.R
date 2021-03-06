`trecase.sex.A` = 
function(yi, ind.lst, X, sex, ni, ni0, xs, start, iphi=1, theta=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  par0 = start
  repeat{
    triali = triali + 1
    tag = tryCatch({
      iphi_i = iphi
      theta_i = theta
      log.lik0 = ll.jRCI.sex.A(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i)
      for(i in 1:maxiter){
        out = optim(par0, ll.jRCI.sex.A, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i)
        log.lik1 = ll.jRCI.sex.A(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i)    
        older = c(rep(par0[1:3], each=2), par0[(4:5)], 0, par0[6], 0, log.lik1, iphi_i, theta_i)
        par0 = out$par
        iphi_i = optimize(f=ll.tRCI.iphi.A, c(tech.ctrl$iphi_l, tech.ctrl$iphi_u), 
                       yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, betas=older[c(1:8, 10)])$minimum
        theta_i = optimize(f=ll.aRC.theta.A, c(tech.ctrl$theta_l, tech.ctrl$theta_u), 
              bs1=older[c(3, 5)], bs2=older[c(4, 6)], ni=ni, ni0=ni0, twosex=TRUE, sex=sex, xs=xs)$minimum    
        log.lik1 = ll.jRCI.sex.A(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i)    
        newer = c(rep(par0[1:3], each=2), par0[(4:5)], 0, par0[6], 0, log.lik1, iphi_i, theta_i)      
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
