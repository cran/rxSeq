`ll.tRCI.iphi.X` = 
function(par0, yi, ind.lst, X, twosex, betas, l.tau.r, l.tau.a){
  iphi = par0
  return(ll.tRCI.X(par0=betas, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi, l.tau.r=l.tau.r, l.tau.a=l.tau.a))
}
