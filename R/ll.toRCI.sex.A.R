`ll.toRCI.sex.A` = 
function(par0, yi, ind.lst, X, iphi){
  par=c(par0[1], par0)
  return(ll.toRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi))
}
