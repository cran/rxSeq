`trecase.sex.b0.X` = 
function(yi, ind.lst, X, ni, ni0, xs, l.tau.r, l.tau.a, 
   start, iphi=1, theta=1, maxiter=100, eps=1E-3, tech.ctrl){
  triali = 0
  par0 = start
  repeat{
    triali = triali + 1
    tag = tryCatch({
      iphi_i = iphi
      theta_i = theta
      log.lik0 = ll.jRCI.sex.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
      for(i in 1:maxiter){
        out = optim(par0, ll.jRCI.sex.b0.X, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)
        log.lik1 = ll.jRCI.sex.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)    
        older = c(par0[1], par0, log.lik1, iphi_i, theta_i)   
        par0 = out$par
        iphi_i = optimize(f=ll.tRCI.iphi.X, c(tech.ctrl$iphi_l, tech.ctrl$iphi_u), 
                       yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, betas=older[1:9], l.tau.r=l.tau.r, l.tau.a=l.tau.a)$minimum
        theta_i = optimize(f=ll.aRC.theta.X, c(tech.ctrl$theta_l, tech.ctrl$theta_u), 
              bs=older[c(3, 4)], ni=ni, ni0=ni0, xs=xs, l.tau.r=l.tau.r)$minimum    
        log.lik1 = ll.jRCI.sex.b0.X(par0, yi=yi, ind.lst=ind.lst, X=X, 
            ni=ni, ni0=ni0, xs=xs, iphi=iphi_i, theta=theta_i, l.tau.r=l.tau.r, l.tau.a=l.tau.a)    
        newer = c(par0[1], par0, log.lik1, iphi_i, theta_i)   
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
