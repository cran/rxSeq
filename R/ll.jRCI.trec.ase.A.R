`ll.jRCI.trec.ase.A` = 
function(par, yi, ind.lst, X, twosex, sex=NULL, ni, ni0, xs, iphi, theta)
{
  if(twosex){
    seq = 1:2
    logL = ll.aRC(par0=c(par[1], par[3]), ni[sex == 1], ni0[sex == 1], xs[sex == 1], theta)
    logL = logL + ll.aRC(par0=c(par[2], par[4]), ni[sex == -1], ni0[sex == -1], xs[sex == -1], theta)  
  }else{
    seq = 1
    logL = ll.aRC(par0=par[1:2], ni, ni0, xs, theta)    
  }
  logL = logL + ll.tRCI.A(par0=c(par[seq], par), yi=yi, ind.lst=ind.lst, X=X, twosex=twosex, iphi=iphi)

  return(logL)
}
                     