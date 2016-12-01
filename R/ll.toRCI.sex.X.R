`ll.toRCI.sex.X` = 
function(par0, yi, ind.lst, X, iphi, l.tau.r, l.tau.a){
  par=c(par0[1], par0)
  return(ll.toRCI.X(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a))
}
