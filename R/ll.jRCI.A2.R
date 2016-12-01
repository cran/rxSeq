`ll.jRCI.A2` = 
function(par, yi, ind.lst, X, twosex, sex, ni, ni0, xs, iphi, theta){
  logL = ll.tRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)
  logL = logL + ll.aRC(par0=par[c(3, 5)], ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
  logL = logL + ll.aRC(par0=par[c(4, 6)], ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)

  return(logL)
}
