`ll.aRC.theta.A` = 
function(par0, bs1, bs2=NULL, ni, ni0, twosex, sex, xs){
  theta = par0
  if(twosex){
    return(ll.aRC(par0=bs1, ni=ni[sex == 1], ni0=ni0[sex == 1], xs=xs[sex == 1], theta=theta) + 
           ll.aRC(par0=bs2, ni=ni[sex ==  -1], ni0=ni0[sex ==  -1], xs=xs[sex ==  -1], theta=theta))
  }else{
    return(ll.aRC(par0=bs1, ni=ni, ni0=ni0, xs=xs, theta=theta))
  }
}