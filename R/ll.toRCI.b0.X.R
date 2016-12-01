`ll.toRCI.b0.X` = 
function(par0, yi, ind.lst, X, twosex, iphi, l.tau.r, l.tau.a){
  if(twosex){
    par=c(0, 0, par0)
  }else{
    par=c(0, par0)
  }

  return(ll.toRCI.X(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a))
}

