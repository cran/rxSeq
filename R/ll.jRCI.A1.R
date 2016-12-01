`ll.jRCI.A1` = 
function(par, yi, ind.lst, X, twosex, ni, ni0, xs, iphi, theta){
  logL = ll.tRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)
  logL = logL + ll.aRC(par0=par[c(2, 3)], ni, ni0, xs, theta)  

  return(logL)
}
