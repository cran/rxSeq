`ll.toRCI.A` = 
function(par0, yi, ind.lst, X, twosex, iphi){
  if(twosex){  
    b0f   = par0[1]
    b0m   = par0[2]
    betas = par0[-(1:2)]
  }else{
    b0f = par0[1]
    betas = par0[-1]
  }  
  etas = rep(0, length(yi)) #base will be eta3 = 0 - AxA, fem and Mal
  ends =  - log(2)
  etas[ind.lst[[4]]] = b0f #BxB, fem
  etas[ind.lst[[1]]] = log1p(exp(b0f)) + ends#AxB, fem
  etas[ind.lst[[2]]] = log1p(exp(b0f)) + ends#BxA, fem
  if(twosex){
    etas[ind.lst[[5]]] = log1p(exp(b0m)) + ends#AxB, mal
    etas[ind.lst[[6]]] = log1p(exp(b0m)) + ends#BxA, mal
    etas[ind.lst[[8]]] = b0m     #BxB, mal
  }
  lmu = c(X %*% betas) + etas

  logL = - loglikNB(iphi=iphi, lmu=lmu, y=yi)
  return(logL)
}
