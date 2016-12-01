`ll.tRCI.A` = 
function(par0, yi, ind.lst, X, twosex, iphi){
  if(twosex){  
    b0f = par0[1]
    b1f = par0[5]
    b0m = par0[2]
    b1m = par0[6]
    betas = par0[ - (1:6)]
  }else{
    b0f = par0[1]
    b1f = par0[3]
    betas = par0[ - (1:3)]
  }  
  etas = rep(0, length(yi)) #base will be eta3 = 0 - AxA, fem
  ends = -log1p(exp(-b1f))
  etas[ind.lst[[4]]] =  b0f #BxB, fem
  etas[ind.lst[[1]]] = -b1f  + log1p(exp(b0f + b1f)) + ends#AxB, fem
  etas[ind.lst[[2]]] =  log1p(exp(b0f - b1f)) + ends#BxA, fem
  if(twosex){
    etas[ind.lst[[5]]] = -b1m + log1p(exp(b0m + b1m)) + ends#AxB, mal
    etas[ind.lst[[6]]] =        log1p(exp(b0m - b1m)) + ends#BxA, mal
    etas[ind.lst[[7]]] =        log1p(exp(-b1m))      + ends#AxA, mal
    etas[ind.lst[[8]]] =  b0m + log1p(exp(-b1m))      + ends#BxB, mal
  }
  lmu = c(X %*% betas) + etas

  logL = -loglikNB(iphi=iphi, lmu=lmu, y=yi)
  return(logL)
}