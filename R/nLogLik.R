`nLogLik` = 
function(res, rc, genei, hessian=FALSE){
  np = ncol(res$coef.full)
  if(rc$model=="full"){
    coef = res$coef.full[genei,1:(np-3)]
    phi = res$coef.full[genei,np-1]
    theta = res$coef.full[genei,np]
    if(rc$chrom=="auto"){
      nll = nLogLik.trecase.A(coef=coef, theta=theta, phi=phi, rc=rc, genei=genei, hessian=hessian)
    }else{
      nll = nLogLik.trecase.X(coef=coef, theta=theta, phi=phi, rc=rc, genei=genei, hessian=hessian)
    }
  }else{
    coef = res$coef.full[genei,1:(np-2)]
    phi = res$coef.full[genei,np]
    if(rc$chrom=="auto"){
      nll = nLogLik.trec.A(coef=coef, phi=phi, rc=rc, genei=genei, hessian=hessian)
    }else{
      nll = nLogLik.trec.X(coef=coef, phi=phi, rc, genei=genei, hessian=hessian)
    }
  }
    
  return(nll)
}

