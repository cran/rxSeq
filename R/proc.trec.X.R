`proc.trec.X` = 
function(rc, hessian=FALSE){
  index = rc$index 
  y = rc$y
  kappas = rc$kappas
  tausB = rc$tausB
  genes.switch = rc$genes.switch
  geneids = rc$geneids

  if(is.null(geneids)){geneids = 1:nrow(y)}
  if(nrow(y) != length(geneids)){
    stop("geneids length is not equal to the number of rows of y")
  }
  if(ncol(y) != length(index)){
    stop("number of columns of y should match length of index")
  }  
  na = sum(index %in% 1:2)
  #genes.switch="ENSMUSG00000086503"

  twosex = max(index)>4
  if(twosex){
    if(!all((1:8) %in% index)){
      stop("package expects to have both sexes of all 4 groups: AB, BA, AA, BB")
    }
    ind.lst = lapply(1:8, function(i){which(index %in% i)})
  }else{
    if(!all((1:4) %in% index)){
      stop("package expects to have all 4 groups: AB, BA, AA, BB")
    }
    ind.lst = lapply(1:4, function(i){which(index %in% i)})
  }  
  if(is.null(kappas)){
    kappas = colSums(y)
  }else{
    if(ncol(y) != length(kappas)){
      stop("number of columns of y should match length of kappas")
    }
  }
  if(is.null(tausB)){
    F1.f = which(index %in% c(1, 2))
    tausB = rep(.5,length(F1.f))
    l.tau.b = l.tau.a = rep(log(.5), length(F1.f))
    l.tau.r = rep(0, length(F1.f))
    warning("Xce effect was not provided, program assumes there is no Xce effect.")
  }else{
    if(na != length(tausB)){
      stop("number of columns of n should match length of tausB")
    }
    l.tau.a = log(1 - tausB); l.tau.b = log(tausB)
    l.tau.r = log(tausB) - l.tau.a      
  }

  l.tau.a = matrix(c(l.tau.a, l.tau.b), ncol=2)
  l.tau.r = matrix(c(l.tau.r, -l.tau.r), ncol=2)


  #maxtrial = 10;
  eps.sh = 1E-3;maxiter.sh = 100
  eps = 1E-3;maxiter = 100
  iphi = 2;theta = .5

  ngenes = length(geneids)
  #genes.switch = c("ENSMUSG00000087640", "ENSMUSG00000086503", "ENSMUSG00000079481", "ENSMUSG00000079444")
  ind = rep(1, ngenes)
  if(length(genes.switch)>0){
    switch.ind = unique(unlist(sapply(genes.switch, grep, geneids)))
    ind[switch.ind] = 2
    message("found ", length(switch.ind), " genes to switch Xce effect:")
    for(i in switch.ind){
      message(geneids[i])
    }
  }
  
  nsamples = length(index)

  intercept = poo = dev = dom = rep(0, length(index))
  intercept = intercept + 1
  dom[index %in% c(1, 2)] = 1
  dev[index %in% c(5, 6)] = 1  
  poo[index == 1] = 1
  poo[index == 2] = -1
  if(twosex){
    sex.tot = intercept;sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dev, poo)
    X.sex = X[, -3]  
    X.dom = X[, -4]
    X.dev = X[, -5]
    X.poo = X[, -6]

    full.model = 2 + ncol(X)
    #add.df = 3; poo.df = 1; dom.df = 1; add.trc.ase = 1; ase.add = 1; sex.df = 2; sex.add.df = 1; dev.df = 1
    df=c(2, 1, 1, 2, 1, 1)
    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom", "pval_sex", "pval_sex.add", "pval_dev")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 2)
    colnames(coef.full) = c("add_F_trc", "add_M_trc", "interc", "kappa", "sex", "dom", "dev.dom", "poo", "ll", "phi")
  }else{
    X = cbind(intercept, kappas, dom, poo)
    X.dom = X[, -3]
    X.poo = X[, -4]

    full.model = 1 + ncol(X)
    #add.df = 3; poo.df = 1; dom.df = 1; add.trc.ase = 1
    df=c(1, 1, 1)

    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 2)
    colnames(coef.full) = c("add_trc", "interc", "kappa", "dom", "poo", "ll", "phi")
  }
  ll.ind = full.model + 1;phi.ind = ll.ind + 1
  
  rownames(pvals) = geneids
  rownames(coef.full) = geneids  

  coef.add =  coef.poo = coef.dom = coef.full
  coef.sex = coef.sex.add = coef.dev.dom = coef.add
  message("processing ", ngenes, " genes")
  errorlist = NULL
  j = 1
  if(twosex){
    for(j in 1:ngenes){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #6. sex effect: b0F' - b0M'=b0F - b0M=b1F - b1M=beta2=beta4=0
      inits.sex = c(0, -5, 1, 0, 0, 0)
      res.sex = trec.sex.X(yi=y[j, ], ind.lst=ind.lst, X=X.sex, 
               l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
               start=inits.sex, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
  
      if(is.null(res.sex)){
        inits.poo = c(0, 0, -5, 1, 0, 0, 0)
        inits.add = c(-5, 1, 0, 0, 0, 0)
        inits.dom = c(0, 0, -5, 1, 0, 0, 0)
        inits.dev = c(0, 0, -5, 1, 0, 0, 0)
        inits.sex.add = c(0, -5, 1, 0, 0, 0)
      }else{
        inits.add = res.sex[3:8]
        inits.poo = res.sex[1:7]
        inits.dom = res.sex[c(1:5, 7:8)]
        inits.dev = res.sex[c(1:6, 8)]
        inits.sex.add = res.sex[2:8]
      }
      #2. POE
      res.poo = trec.b1.X(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
               l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
               start=inits.poo, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #1. Strain
      res.add = trec.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
              l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
              start=inits.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #3. Dominance
      res.dom = trec.dom.X(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
              l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
              start=inits.dom, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #8. Deviating dominance
      res.dev = trec.dev.X(yi=y[j, ], ind.lst=ind.lst, X=X.dev, 
              l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
              start=inits.dev, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   

      #7. sex specific addain: b0F' - b0M'=b0F - b0M=0
      res.sex.add = trec.sex.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, 
               l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
               start=inits.sex.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      
      #find better start for full    
      #"pval_add", "pval_poo", "pval_dom", "pval_sex", "pval_sex.add", "pval_dev"
      lls = list(res.add, res.poo, res.dom, res.sex, res.sex.add, res.dev)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)                         
      #fit full
      res.full = trec.full.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      #refit add
      inits.add = res.full[3:8]
      inits.poo = res.full[1:7]
      inits.dom = res.full[c(1:5, 7:8)]
      inits.dev = res.full[c(1:6, 8)]
      inits.sex = c(mean(res.full[1:2]), res.full[c(3:4, 6:8)])
      inits.sex.add = c(mean(res.full[1:2]), res.add[3:8])

      res.add1 = trec.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      res.poo1 = trec.b1.X(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.poo, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      res.dom1 = trec.dom.X(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.dom, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit sex
      res.sex1 = trec.sex.X(yi=y[j, ], ind.lst=ind.lst, X=X.sex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.sex, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #5. refit sex specific addain
      res.sex.add1 = trec.sex.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.sex.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #7. sex specific Dominance
      res.dev1 = trec.dev.X(yi=y[j, ], ind.lst=ind.lst, X=X.dev, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.dev, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind, trecase=FALSE)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind, trecase=FALSE)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind, trecase=FALSE)
      res.sex = get.better(lst=list(res.sex, res.sex1), ind=ll.ind, trecase=FALSE)
      res.sex.add = get.better(lst=list(res.sex.add, res.sex.add1), ind=ll.ind, trecase=FALSE)
      res.dev = get.better(lst=list(res.dev, res.dev1), ind=ll.ind, trecase=FALSE)
  
      #find better start for full
                         #c("pval_add", "pval_poo", "pval_dom", "pval_sex", "pval_sex.add", "pval_dev")
      lls = list(res.full, res.add, res.poo, res.dom, res.sex, res.sex.add, res.dev)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)
      #fit full
      res.full1 = trec.full.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind, trecase=FALSE)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
      if(!is.null(res.sex)){coef.sex[j, ] = res.sex}
      if(!is.null(res.sex.add)){coef.sex.add[j, ] = res.sex.add}
      if(!is.null(res.dev)){coef.dev.dom[j, ] = res.dev}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    coef.sex[, phi.ind] = 1/coef.sex[, phi.ind]
    coef.sex.add[, phi.ind] = 1/coef.sex.add[, phi.ind]
    coef.dev.dom[, phi.ind] = 1/coef.dev.dom[, phi.ind]
    #c("pval_add", "pval_poo", "pval_dom", "pval_sex", "pval_sex.add", "pval_dev.dom")
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, 
    coef.dom=coef.dom, coef.sex=coef.sex, coef.sex.add=coef.sex.add, coef.dev.dom=coef.dev.dom, 
    errorlist=errorlist)
  }else{
    for(j in 1:ngenes){
    #for(j in 1:2){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #1. Strain
      inits.add = c(-5, 1, 0, 0)
      #yi=y[j, ];sex=sex.ase; ni=nis[j, ];ni0=ni0B[j, ];start=inits.add
      res.add = trec.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
              l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
              start=inits.add, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)     
      if(is.null(res.add)){
        inits.poo = c(0, -5, 1, 0)
        inits.dom = c(0, -5, 1, 0)
      }else{
        inits.poo = res.add[1:4]
        inits.dom = res.add[c(1:3, 5)]
      }
      #2. POE
      res.poo = trec.b1.X(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
               l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
               start=inits.poo, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #3. Dominance
      res.dom = trec.dom.X(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
              l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
              start=inits.dom, iphi=iphi, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      
      #find better start for full    
      lls = list(res.add, res.poo, res.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)                         
      #fit full
      res.full = trec.full.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      #refit add
      inits.add = res.full[2:5]
      inits.poo = res.full[c(1:2, 4:5)]
      inits.dom = res.full[1:4]

      res.add1 = trec.b0.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.add, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      res.poo1 = trec.b1.X(yi=y[j, ], ind.lst=ind.lst, X=X.poo, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.poo, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      res.dom1 = trec.dom.X(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=inits.dom, iphi=iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind, trecase=FALSE)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind, trecase=FALSE)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind, trecase=FALSE)
  
      #find better start for full
      lls = list(res.full, res.add, res.poo, res.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind, trecase=FALSE)
      #fit full
      res.full1 = trec.full.X(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, 
             l.tau.r=l.tau.r[, ind[j]], l.tau.a=l.tau.a[, ind[j]], 
             start=new.inits$inits, iphi=new.inits$iphi_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind, trecase=FALSE)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, 
                   coef.dom=coef.dom, errorlist=errorlist)  
  }  
  if(hessian){
    if(ngenes>100)warning("This option is recommended only for the extra analysis on a subset of genes.")    
    hess.lst = as.list(rep(NA, ngenes))
    for(j in 1:ngenes){
      tag = tryCatch({    
        nll = nLogLik.trec.X(coef=coef.full[j, 1:(ll.ind-1)], phi=coef.full[j, phi.ind], 
          rc, j, hessian=TRUE)
        hess.lst[[j]] = nll$hess
        0
      }, error=function(e) {
        warning("a problem with hessian")
        1
      })
      names(hess.lst) = geneids
    }
    res$hess.lst = hess.lst
  }  
  return(res)
}

