`ll.tRCI.iphi.A` = 
function(par0, yi, ind.lst, X, twosex, betas){
  iphi = par0
  return(ll.tRCI.A(par0=betas, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi))
}