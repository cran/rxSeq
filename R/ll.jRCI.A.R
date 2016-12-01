`ll.jRCI.A` = 
function(par, yi, ind.lst, X, twosex, sex, ni, ni0, xs, iphi, theta){
  logL = ll.tRCI.A(par0=par, yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)
  if(twosex){
    ase = c(3, 5)
    logL = logL + ll.aRC(par0=par[ase], ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=par[ase + 1], ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  }else{
    ase = c(2, 3)
    logL = logL + ll.aRC(par0=par[ase], ni, ni0, xs, theta)
  }

  return(logL)
}
