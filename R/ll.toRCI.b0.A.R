`ll.toRCI.b0.A` = 
function(par0, yi, ind.lst, X, twosex, iphi){
  if(twosex){
    par=c(0, 0, par0)
  }else{
    par=c(0, par0)
  }

  return(ll.toRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi))
}
