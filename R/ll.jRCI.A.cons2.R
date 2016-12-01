`ll.jRCI.A.cons2` = 
function(par, yi, ind.lst, X, sex, ni, ni0, xs, iphi, theta){
  par1 = c(par[1:2], par)
  logL = ll.tRCI.A(par0=par1, yi=yi, ind.lst=ind.lst, X=X, twosex=TRUE, iphi=iphi)
  logL = logL + ll.aRC(par0=par1[c(3, 5)], ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
  logL = logL + ll.aRC(par0=par1[c(4, 6)], ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)

  return(logL)
}