`proc.trecase.A` = 
function(rc, hessian=FALSE){
  index = rc$index 
  y = rc$y
  n = rc$n 
  n0B = rc$n0B
  kappas = rc$kappas
  geneids = rc$geneids
  
  if(nrow(y) != nrow(n)){
    stop("number of rows of n should match number of rows of y")
  }
  if(nrow(n0B) != nrow(n)){
    stop("number of rows of n should match number of rows of y0B")
  }
  if(is.null(geneids)){geneids = 1:nrow(y)}
  if(nrow(y) != length(geneids)){
    stop("geneids length mismatches with number of rows of n")
  }
  if(ncol(y) != length(index)){
    stop("number of rows of y should match length of index")
  }  

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

  #maxtrial = 10;
  eps.sh = 1E-2; maxiter.sh = 20
  eps = 1E-3; maxiter = 100
  iphi = 2; theta = .5

  ngenes = length(geneids)
  nsamples = length(index)
  intercept = rep(1, nsamples)
  dom = as.numeric(index %in% c(1, 2, 5, 6))
  xs = index[1:ncol(n)]; xs[xs == 5] = 1;xs[xs %in% c(2, 6)] = -1
  if(twosex){
    sex.tot = intercept;sex.tot[index>4] = -1
    X = cbind(intercept, kappas, sex.tot, dom, dom*sex.tot)
    X.sex.dom = X[, -5]
    X.dom = X.sex.dom[, -4]
    X.sex = X[, c(-3, -5)]  

    sex.ase = sex.tot[1:ncol(n)]
    full.model = 6 + ncol(X)
    df=c(4, 2, 2, 2, 2, 5, 2, 1, 1)
    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom", "pval_same", "pval_ase.add", "pval_sex", "pval_sex.add", "pval_sex.poo", "pval_sex.dom")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 3)
    colnames(coef.full) = c("add_F_trc", "add_M_trc", "add_F_ase", "add_M_ase", "poo_F", "poo_M", "interc", "kappa", "sex", "dom", "sex.dom", "ll", "phi", "theta")
  }else{
    X = cbind(intercept, kappas, dom)
    X.dom = X[, -3]

    sex.ase = NULL
  
    full.model = 3 + ncol(X)
    df=c(2, 1, 1, 1, 1)
    pvals = matrix(NA, nrow=ngenes, ncol=length(df))
    colnames(pvals) = c("pval_add", "pval_poo", "pval_dom", "pval_same", "pval_ase.add")
    coef.full = matrix(NA, nrow=ngenes, ncol=full.model + 3)
    colnames(coef.full) = c("add_trc", "add_ase", "poo", "interc", "kappa", "dom", "ll", "phi", "theta")
  }
  ll.ind = full.model + 1; phi.ind = ll.ind + 1; theta.ind = ll.ind + 2
  rownames(pvals) = geneids
  rownames(coef.full) = geneids  

  coef.ase.add = coef.add =  coef.poo = coef.dom = coef.full
  coef.sex = coef.sex.add = coef.sex.poo = coef.sex.dom = coef.same = coef.add
  message("processing ", ngenes, " genes")
  errorlist = NULL
  if(twosex){
    for(j in 1:ngenes){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #4. sex effect: b0F' - b0M'=b0F - b0M=b1F - b1M=beta2=beta4=0
      #yi=y[j, ];X=X.sex;sex=sex.ase;ni=n[j, ];ni0=n0B[j, ];start=inits.sex;maxiter=maxiter.sh;eps=eps.sh
      inits.sex = c(0, 0, 0, -5, 1, 0)
      res.sex = trecase.sex.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex, sex=sex.ase, 
               ni=n[j, ], ni0=n0B[j, ], xs=xs, 
               start=inits.sex, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
  
      if(is.null(res.sex)){
        inits.poo = c(0, 0, 0, 0, -5, 1, 0, 0, 0)
        inits.add = c(0, 0, -5, 1, 0, 0, 0)
        inits.dom = c(0, 0, 0, 0, 0, 0, -5, 1, 0)
      }else{
        inits.add = res.sex[5:11]
        inits.poo = res.sex[c(1:4, 7:11)]
        inits.dom = res.sex[1:9]
      }
      #2. POE
      res.poo = trecase.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
               ni=n[j, ], ni0=n0B[j, ], xs=xs, 
               start=inits.poo, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #1. Strain
      res.add = trecase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
              ni=n[j, ], ni0=n0B[j, ], xs=xs, 
              start=inits.add, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #3. Dominance
      res.dom = trecase.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, sex=sex.ase, 
              ni=n[j, ], ni0=n0B[j, ], xs=xs, 
              start=inits.dom, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #5. sex specific addain: b0F' - b0M'=b0F - b0M=0
      if(is.null(res.add)){inits.sex.add = c(0, 0, 0, 0, -5, 1, 0, 0, 0)}else{inits.sex.add = res.add[3:11]}
      res.sex.add = trecase.sex.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, sex=sex.ase, 
               ni=n[j, ], ni0=n0B[j, ], xs=xs, 
               start=inits.sex.add, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #6. sex specific poo: b1F=b1M
      if(is.null(res.poo)){inits.sex.poo = c(0, 0, 0, 0, 0, -5, 1, 0, 0, 0)}else{inits.sex.poo = res.poo[c(1:5, 7:11)]}
      res.sex.poo = trecase.sex.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, sex=sex.ase, 
               ni=n[j, ], ni0=n0B[j, ], xs=xs, 
               start=inits.sex.poo, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #7. sex specific Dominance
      if(is.null(res.dom)){inits.sex.dom = c(0, 0, 0, 0, 0, 0, -5, 1, 0, 0)}else{inits.sex.dom = res.dom[1:10]}
      res.sex.dom = trecase.sex.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.dom, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.sex.dom, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #8. trec and ace
      if(is.null(res.sex.dom)){inits.same = c(0, 0, 0, 0, -5, 1, 0, 0, 0)}else{
        inits.same = c(mean(res.sex.dom[c(1, 3)]), mean(res.sex.dom[c(2, 4)]), res.sex.dom[5:11])
      }
      res.same = trecase.trec.ase.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.same, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #9. ase only b0
      if(is.null(res.sex.add)){inits.ase.b0 = c(0, 0, 0, 0, -5, 1, 0, 0, 0)}else{inits.ase.b0 = res.sex.add[c(1:2, 5:11)]}
      res.ase.add = trecase.ase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
            ni=n[j, ], ni0=n0B[j, ], xs=xs, 
            start=inits.ase.b0, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      
      #find better start for full    
      #c("pval_add", "pval_poo", "pval_dom", "pval_same", "pval_ase.add", "pval_sex", "pval_sex.add", "pval_sex.poo", "pval_sex.dom")      
      lls = list(res.add, res.poo, res.dom, res.same, res.ase.add, res.sex, res.sex.add, res.sex.poo, res.sex.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind)                         
      #fit full
      res.full = trecase.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=new.inits$inits, iphi=new.inits$iphi_i, theta=new.inits$theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      theta_i = res.full[theta.ind]
      #refit add
      inits.add = res.full[5:11]
      inits.poo = res.full[c(1:4, 7:11)]
      inits.dom = res.full[1:9]
      inits.sex = c(mean(res.full[1:2]), mean(res.full[3:4]), mean(res.full[5:6]), res.full[c(7:8, 10)])
      inits.sex.add = c(mean(res.full[1:2]), mean(res.full[3:4]), res.add[5:11])
      inits.sex.poo = c(res.full[1:4], mean(res.full[5:6]), res.add[7:11])
      inits.sex.dom = res.full[1:10]
      inits.same = c(mean(res.full[c(1, 3)]), mean(res.full[c(2, 4)]), res.full[5:11])
      inits.ase.add = res.full[c(1:2, 5:11)]
      res.add1 = trecase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.add, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      res.poo1 = trecase.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.poo, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      res.dom1 = trecase.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.dom, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit sex
      res.sex1 = trecase.sex.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.sex, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #5. refit sex specific addain
      res.sex.add1 = trecase.sex.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.sex.add, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #6. sex specific poo: b1F=b1M
      res.sex.poo1 = trecase.sex.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.sex.poo, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #7. sex specific Dominance
      res.sex.dom1 = trecase.sex.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.sex.dom, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.sex.dom, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #8. trec and ace
      res.same1 = trecase.trec.ase.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.same, iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #9. ase only b0 
      res.ase.add1 = trecase.ase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.ase.add, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind)
      res.sex = get.better(lst=list(res.sex, res.sex1), ind=ll.ind)
      res.sex.add = get.better(lst=list(res.sex.add, res.sex.add1), ind=ll.ind)
      res.sex.poo = get.better(lst=list(res.sex.poo, res.sex.poo1), ind=ll.ind)
      res.sex.dom = get.better(lst=list(res.sex.dom, res.sex.dom1), ind=ll.ind)
      res.same = get.better(lst=list(res.same, res.same1), ind=ll.ind)
      res.ase.add = get.better(lst=list(res.ase.add, res.ase.add1), ind=ll.ind)
  
      #find better start for full
      #                 c("pval_add", "pval_poo", "pval_dom", "pval_same", "pval_ase.add", "pval_sex", "pval_sex.add", "pval_sex.poo", "pval_sex.dom")
      lls = list(res.full, res.add, res.poo, res.dom, res.same, res.ase.add, res.sex, res.sex.add, res.sex.poo, res.sex.dom)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind)
      #fit full
      res.full1 = trecase.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=new.inits$inits, iphi=new.inits$iphi_i, theta=new.inits$theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
      if(!is.null(res.sex)){coef.sex[j, ] = res.sex}
      if(!is.null(res.sex.add)){coef.sex.add[j, ] = res.sex.add}
      if(!is.null(res.sex.poo)){coef.sex.poo[j, ] = res.sex.poo}
      if(!is.null(res.sex.dom)){coef.sex.dom[j, ] = res.sex.dom}
      if(!is.null(res.same)){coef.same[j, ] = res.same}
      if(!is.null(res.ase.add)){coef.ase.add[j, ] = res.ase.add}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    coef.same[, phi.ind] = 1/coef.same[, phi.ind]
    coef.ase.add[, phi.ind] = 1/coef.ase.add[, phi.ind]    
    coef.sex[, phi.ind] = 1/coef.sex[, phi.ind]
    coef.sex.add[, phi.ind] = 1/coef.sex.add[, phi.ind]
    coef.sex.poo[, phi.ind] = 1/coef.sex.poo[, phi.ind]
    coef.sex.dom[, phi.ind] = 1/coef.sex.dom[, phi.ind]
    #c("pval_add", "pval_poo", "pval_dom", "pval_same", "pval_ase.add", "pval_sex", "pval_sex.add", "pval_sex.poo", "pval_sex.dom")
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, 
    coef.dom=coef.dom, coef.same=coef.same, coef.ase.add=coef.ase.add, coef.sex=coef.sex, 
    coef.sex.add=coef.sex.add, coef.sex.poo=coef.sex.poo, coef.sex.dom=coef.sex.dom, errorlist=errorlist)
  }else{
    for(j in 1:ngenes){
      if(j%%100 == 0){
        message(j, "th gene")
      }
      #fit the smallest first 
      #1. Strain
      inits.add = c(0, -5, 1, 0)
      #yi=y[j, ];sex=sex.ase; ni=n[j, ];ni0=n0B[j, ];start=inits.add
      res.add = trecase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
              ni=n[j, ], ni0=n0B[j, ], xs=xs, 
              start=inits.add, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)     
      if(is.null(res.add)){
        inits.poo = c(0, 0, -5, 1, 0)
        inits.dom = c(0, 0, 0, -5, 1)
      }else{
        inits.poo = res.add[c(1:2, 4:6)]
        inits.dom = res.add[1:5]
      }
      #POE
      res.poo = trecase.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
               ni=n[j, ], ni0=n0B[j, ], xs=xs, 
               start=inits.poo, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      #Dominance
      res.dom = trecase.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, sex=sex.ase, 
              ni=n[j, ], ni0=n0B[j, ], xs=xs, 
              start=inits.dom, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #trec and ace
      if(is.null(res.poo)){inits.same = c(0, 0, -5, 1, 0)}else{
        inits.same = c(mean(res.poo[1:2]), res.poo[3:6])
      }
      res.same = trecase.trec.ase.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.same, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)   
      #ase only b0
      if(is.null(res.poo)){inits.ase.add = c(0, 0, -5, 1, 0)}else{inits.ase.add = res.poo[c(1, 3:6)]}
      res.ase.add = trecase.ase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
            ni=n[j, ], ni0=n0B[j, ], xs=xs, 
            start=inits.ase.add, iphi=iphi, theta=theta, maxiter=maxiter.sh, eps=eps.sh, tech.ctrl=rc$tech.ctrl)
      
      #find better start for full    
      lls = list(res.add, res.poo, res.dom, res.same, res.ase.add)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind)                         
      #fit full
      res.full = trecase.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=new.inits$inits, iphi=new.inits$iphi_i, theta=new.inits$theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      if(is.null(res.full)){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next          
      }
      iphi_i = res.full[phi.ind]
      theta_i = res.full[theta.ind]
      #refit add
      inits.add = res.full[3:6]
      res.add1 = trecase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.add, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit poo
      inits.poo = res.full[c(1:2, 4:6)]
      res.poo1 = trecase.b1.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.poo, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #refit dom
      inits.dom = res.full[1:5]
      res.dom1 = trecase.dom.A(yi=y[j, ], ind.lst=ind.lst, X=X.dom, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.dom, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #8. trec and ace
      inits.same = c(mean(res.full[c(1, 2)]), res.full[3:6])
      res.same1 = trecase.trec.ase.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.same, iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      #9. ase only b0 
      inits.ase.add = res.full[c(1, 3:6)]
      res.ase.add1 = trecase.ase.b0.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=inits.ase.add, iphi=iphi_i, theta=theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.add = get.better(lst=list(res.add, res.add1), ind=ll.ind)
      res.poo = get.better(lst=list(res.poo, res.poo1), ind=ll.ind)
      res.dom = get.better(lst=list(res.dom, res.dom1), ind=ll.ind)
      res.same = get.better(lst=list(res.same, res.same1), ind=ll.ind)
      res.ase.add = get.better(lst=list(res.ase.add, res.ase.add1), ind=ll.ind)
  
      #find better start for full
      lls = list(res.full, res.add, res.poo, res.dom, res.same, res.ase.add)
      if(all(sapply(lls, is.null))){
        errorlist = c(errorlist, j) 
        warning("no test ", j)
        next   
      }
      new.inits = get.inits(lst=lls, ind=ll.ind)
      #fit full
      res.full1 = trecase.full.A(yi=y[j, ], ind.lst=ind.lst, X=X, twosex=twosex, sex=sex.ase, 
             ni=n[j, ], ni0=n0B[j, ], xs=xs, 
             start=new.inits$inits, iphi=new.inits$iphi_i, theta=new.inits$theta_i, maxiter=maxiter, eps=eps, tech.ctrl=rc$tech.ctrl)
      res.full = get.better(lst=list(res.full, res.full1), ind=ll.ind)
      coef.full[j, ] = res.full
      tests = test(full=res.full, short=lls[-1], df=df, ind=ll.ind, twosex=twosex, genei=j)
      pvals[j, ] = tests$pvals
      errorlist = c(errorlist, tests$error)
      if(!is.null(res.add)){coef.add[j, ] = res.add}
      if(!is.null(res.poo)){coef.poo[j, ] = res.poo}
      if(!is.null(res.dom)){coef.dom[j, ] = res.dom}
      if(!is.null(res.same)){coef.same[j, ] = res.same}
      if(!is.null(res.ase.add)){coef.ase.add[j, ] = res.ase.add}
    }
    coef.full[, phi.ind] = 1/coef.full[, phi.ind]
    coef.add[, phi.ind] = 1/coef.add[, phi.ind]
    coef.poo[, phi.ind] = 1/coef.poo[, phi.ind]    
    coef.dom[, phi.ind] = 1/coef.dom[, phi.ind]    
    coef.same[, phi.ind] = 1/coef.same[, phi.ind]
    coef.ase.add[, phi.ind] = 1/coef.ase.add[, phi.ind]    
    res = list(pvals=pvals, coef.full=coef.full, coef.add=coef.add, coef.poo=coef.poo, 
    coef.dom=coef.dom, coef.same=coef.same, coef.ase.add=coef.ase.add, errorlist=errorlist)  
  }  
  if(hessian){
    if(ngenes>100)warning("This option is recommended only for the extra analysis on a subset of genes.")    
    hess.lst = as.list(rep(NA, ngenes))
    for(j in 1:ngenes){
      tag = tryCatch({    
        nll = nLogLik.trecase.A(coef=coef.full[j, 1:(ll.ind-1)], 
          phi=coef.full[j, phi.ind], theta=coef.full[j, theta.ind], 
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

