`process` = 
function(rc, hessian=FALSE){
#  tech.ctrl = rc$tech.ctrl
#  globalVariables("tech.ctrl")
  if(!(is.null(rc$n)|is.null(rc$n0B))){
    if(nrow(rc$y) != nrow(rc$n)){
      stop("number of rows of n should match number of rows of y")
    }
    if(nrow(rc$n0B) != nrow(rc$n)){
      stop("number of rows of n should match number of rows of n0B")
    }
  }
  if(is.null(rc$geneids)){rc$geneids = 1:nrow(rc$y)}
  if(nrow(rc$y) != length(rc$geneids)){
    stop("geneids length is not equal to the number of rows of y")
  }
  if(ncol(rc$y) != length(rc$index)){
    stop("number of columns of y should match length of index")
  }  
  if(rc$model=="full"){
    if(rc$chrom=="auto"){
      res = proc.trecase.A(rc, hessian=hessian)
    }else{
      res = proc.trecase.X(rc, hessian=hessian)
    }
  }else{
    if(rc$chrom=="auto"){
      res = proc.trec.A(rc, hessian=hessian)
    }else{
      res = proc.trec.X(rc, hessian=hessian)
    }
  }
    
  return(res)
}

