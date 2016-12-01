`ll.tRCI.X` = 
function(par0, yi, ind.lst, X, twosex, iphi, l.tau.r, l.tau.a){
  if(twosex){  
    b0f = par0[1]
    b0m = par0[2]
    #b0f ase par0[3]
    b1f = par0[4]
    betas = par0[-(1:4)]
  }else{
    b0f = par0[1]
    #b0f ase par0[2]
    b1f = par0[3]
    betas = par0[-(1:3)]
  }
  
  end = log(2) - log1p(exp(b1f))
  etas = rep(0, length(yi))                                                        #AxA, fem - base
  etas[ind.lst[[4]]] = b0f                                                         #BxB, fem
  etas[ind.lst[[1]]] =       log1p(exp(l.tau.r[ind.lst[[1]]] + b0f + b1f)) + l.tau.a[ind.lst[[1]]] + end#AxB, fem
  etas[ind.lst[[2]]] = b1f + log1p(exp(l.tau.r[ind.lst[[2]]] + b0f - b1f)) + l.tau.a[ind.lst[[2]]] + end#BxA, fem
  if(twosex){
    etas[ind.lst[[5]]] = b1f + end                                                 #AxB, mal
    etas[ind.lst[[6]]] = b0m + b1f + end                                           #BxA, mal
    etas[ind.lst[[7]]] = b1f + end                                                 #AxA, mal
    etas[ind.lst[[8]]] = b0m + b1f + end                                           #BxB, mal
  }  

  lmu = c(X%*%betas) + etas

  logL = -loglikNB(iphi=iphi, lmu=lmu, y=yi)
  return(logL)
}
