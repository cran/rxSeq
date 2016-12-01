`ll.aRC.theta.X` = 
function(par0, bs, ni, ni0, xs, l.tau.r){
  theta = par0
  return(ll.aRC.X(par0=bs, ni=ni, ni0=ni0, xs=xs, theta=theta, l.tau.r=l.tau.r))
}
