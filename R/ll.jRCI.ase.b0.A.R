`ll.jRCI.ase.b0.A` = 
function(par, yi, ind.lst, X, twosex, sex, ni, ni0, xs, iphi, theta)
{ 
  if(twosex){
    seq = 1:2
    zeroes = c(0, 0)
    logL = ll.aRC(par0=c(0, par[3]), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=c(0, par[4]), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)
  }else{
    seq = 1
    zeroes = 0
    logL = ll.aRC(par0=c(0, par[2]), ni, ni0, xs, theta)  
  }
  logL = logL + ll.tRCI.A(par0=c(par[seq], zeroes, par[-seq]), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)

  return(logL)
}
