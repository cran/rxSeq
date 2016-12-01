`ll.aRC.theta.A1` = 
function(par0, bs1, bs2=NULL, ni, ni0, xs){
  theta = par0
  return(ll.aRC(par0=bs1, ni=ni, ni0=ni0, xs=xs, theta=theta))
}