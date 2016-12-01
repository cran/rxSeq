`ll.jRCI.A.cons` = 
function(par, yi, ind.lst, X, twosex, sex, ni, ni0, xs, iphi, theta){
  if(twosex){
    par1 = c(par[1:2], par)
    logL = ll.aRC(par0=par1[c(3, 5)], ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=par1[c(4, 6)], ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  }else{
    par1 = c(par[1], par)
    logL = ll.aRC(par0=par1[c(2, 3)], ni, ni0, xs, theta)  
  }
  logL = logL + ll.tRCI.A(par0=par1, yi=yi, ind.lst=ind.lst, X=X, twosex == twosex, iphi=iphi)
  return(logL)
}