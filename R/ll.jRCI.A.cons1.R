`ll.jRCI.A.cons1` = 
function(par, yi, ind.lst, X, sex, ni, ni0, xs, iphi, theta){
  par1 = c(par[1], par)
  logL = ll.tRCI.A(par0=par1, yi=yi, ind.lst=ind.lst, X=X, twosex=FALSE, iphi=iphi)
  logL = logL + ll.aRC(par0=par1[c(2, 3)], ni, ni0, xs, theta)  

  return(logL)
}