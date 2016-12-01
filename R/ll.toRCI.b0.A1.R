`ll.toRCI.b0.A1` = 
function(par0, yi, ind.lst, X, iphi){
  par=c(0, par0)

  return(ll.toRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=FALSE, iphi=iphi))
}
