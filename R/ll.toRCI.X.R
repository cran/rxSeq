`ll.toRCI.X` = 
function(par0, yi, ind.lst, X, twosex, iphi, l.tau.r, l.tau.a){
  if(twosex){  
    b0f = par0[1]
    b0m = par0[2]
    betas = par0[-(1:2)]
  }else{
    b0f = par0[1]
    betas = par0[-1]
  }  
  etas               = rep(0, length(yi))                                             #AxA, fem AxA, mal, AxB, mal #3, 7, 5
  etas[ind.lst[[4]]] = b0f                                                            #BxB, fem                  #4
  etas[ind.lst[[1]]] = log1p(exp(l.tau.r[ind.lst[[1]]] + b0f)) + l.tau.a[ind.lst[[1]]]#AxB, fem                  #1
  etas[ind.lst[[2]]] = log1p(exp(l.tau.r[ind.lst[[2]]] + b0f)) + l.tau.a[ind.lst[[2]]]#BxA, fem                  #2
  if(twosex){
    etas[ind.lst[[6]]] = b0m                                                          #BxA, mal                  #6
    etas[ind.lst[[8]]] = b0m                                                          #BxB, mal                  #8
  }
  lmu = c(X%*%betas) + etas

  logL = -loglikNB(iphi=iphi, lmu=lmu, y=yi)
  return(logL)
}
