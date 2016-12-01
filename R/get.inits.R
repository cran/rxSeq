`get.inits` = 
function(lst, ind, trecase=TRUE){
  lls = Inf
  theta_i = iphi_i = inits = NULL
  for(item in lst){
    if(!is.null(item)){
      if(lls>item[ind]){
        lls = item[ind]
        inits = item[1:(ind-1)]
        iphi_i = item[ind + 1]
        if(trecase)theta_i = item[ind + 2]
      }
    }
  }
  if(lls<Inf){
    if(trecase){
      return(list(inits=inits, iphi_i=iphi_i, theta_i=theta_i, ll=lls))  
    }else{
      return(list(inits=inits, iphi_i=iphi_i, ll=lls))  
    }
  }else{
    return(NULL)
  }
}

