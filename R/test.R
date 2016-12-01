`test` = 
function(full, short, df, ind, twosex, genei=""){
  n = length(short)
  if(twosex){
    hyp=c("str", "poo", "dom", "same", "str.ase", "sex", "sex.str", "sex.poo", "sex.dom")
  }else{
    hyp=c("str", "poo", "dom", "same", "str.ase")
  }
  pvals = rep(NA, n)
  error = NULL
  for(i in 1:n){
    if(is.null(full) | is.null(short[[i]])){
      error = c(error, paste(hyp[i], genei, sep=":"))
    }else{
      lrt = 2 * ((short[[i]])[ind] - full[ind])
      pvals[i] = pchisq(lrt, df=df[i], lower.tail=F)
    }
  }
  return(list(pvals=pvals, error=error))  
}
