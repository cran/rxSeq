`ll.aRC.theta.A2` = 
function(par0, bs1, bs2, ni, ni0, sex, xs){
  theta = par0
  return(ll.aRC(par0=bs1, ni=ni[sex == 1], ni0=ni0[sex == 1], xs=xs[sex == 1], theta=theta) + 
         ll.aRC(par0=bs2, ni=ni[sex == -1], ni0=ni0[sex == -1], xs=xs[sex == -1], theta=theta))
}