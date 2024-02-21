
simpute.als.fit_splr <-
function (y, X=NULL, H=NULL, J = 2, thresh = 1e-05, lambda=0, 
          maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE,
          patience=3, svdH=NULL, return_obj=FALSE, init="naive") {

  if(!inherits(y,"dgCMatrix")) y=as(y,"dgCMatrix")
  irow=y@i
  pcol=y@p
  n <- dim(y)
  m <- n[2]
  n <- n[1]
  nz=nnzero(y)
  initialize_beta = FALSE
  #-------------------------------
  if(is.null(svdH)){
    stopifnot(!is.null(X))
    svdH = reduced_hat_decomp(X, 1e-2)
    J_H = svdH$rank
    svdH = svdH$svdH
    if(trace.it) print(paste("Rank of H is ", J_H))
  }
  #---------------------------------------------------
  warm=FALSE
  clean.warm.start(warm.start)
  if(!is.null(warm.start)){
    #must have u,d and v components
    if(!all(match(c("u","d","v","xbeta.obs"),names(warm.start),0)>0))
      stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=sum(D>0)
    #J = JD
    xbeta.obs = warm.start$xbeta.obs
    if(JD >= J){
      U=warm.start$u[,seq(J),drop=FALSE]
      V=warm.start$v[,seq(J),drop=FALSE]
      Dsq=D[seq(J)]
    }else{
      Dsq=c(D,rep(D[JD],J-JD))
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=fast.svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
    if(any(is.na(warm.start$xbeta.obs))){
      initialize_beta = TRUE
      #xbeta.obs = suvC(svdH$u, t(as.matrix(svdH$v %*% y)),irow,pcol)
    }else{
      xbeta.obs = warm.start$xbeta.obs
    }
  }else{
    stopifnot(init %in% c("random", "naive"))
    if(init == "random"){
      
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=fast.svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VDsq and AB'=U Dsq V^T
      xbeta.obs = suvC(svdH$u, t(as.matrix(svdH$v %*% y)),irow,pcol)
    }else if(init == "naive"){
      Y_naive = as.matrix(y)
      yobs = ! is.na(Y_naive)
      Y_naive = naive_MC(Y_naive)
      naive_fit <-  svdH$u %*% (svdH$v  %*% Y_naive)
      xbeta.obs <- naive_fit[yobs]
      naive_fit <- Y_naive - naive_fit
      naive_fit <- propack.svd(as.matrix(naive_fit), J)
      U = naive_fit$u
      V = naive_fit$v
      Dsq = naive_fit$d
    }
  }
  
  #----------------------------------------
  S=y
  HU = svdH$u %*% (svdH$v %*% U)
  ratio <- 1
  iter <- 0
  counter = 0
  best_score = Inf
  best_iter = NA
  if(return_obj) obj.l <- rep(NA, maxit)
  #-----------------------------------------------------------
  if(initialize_beta){
    VDsq=UD(V,Dsq,m)
    S@x = y@x - suvC(U,VDsq,irow,pcol)
    part1 = suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol)
    part2 = suvC(HU,VDsq, irow, pcol)
    xbeta.obs = part1 + part2 + suvC(svdH$u, t(as.matrix(svdH$v %*% y)),irow,pcol)
  }
  #----------------------------------------
  while ((ratio > thresh)&(iter < maxit)) {
    iter <- iter + 1
    U.old=U
    V.old= V
    Dsq.old=Dsq
    #---------------------------------
    ## U step # S is yplus
    if(iter>1|warm){

      VDsq=UD(V,Dsq,m)
      M_obs = suvC(U,VDsq,irow,pcol)
      S@x = y@x - M_obs - xbeta.obs  
  
    }else {
      VDsq=matrix(0,m,J) 
      S@x = y@x - xbeta.obs
    }
    
    IHU =  U - HU
    # Compute beta estimates
    part1 = suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol)
    part2 = suvC(HU,VDsq, irow, pcol)
    xbeta.obs = part1 + part2 + xbeta.obs
    # update B
    B = as.matrix(t(S) %*% IHU) + (VDsq %*% (t(U) %*% IHU))  
    if(lambda>0) B = t(t(B) *(Dsq/(Dsq+lambda))) 
    Bsvd=fast.svd(B)
    V = Bsvd$u      
    Dsq = Bsvd$d
    U = U%*% (Bsvd$v)
    #------------------------------------
    # V step
    ## Compute beta estimates again
    HU = svdH$u %*% (svdH$v %*% U)
    part1 = suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol)
    part2 = suvC(HU,VDsq, irow, pcol)
    xbeta.obs = part1 + part2 + xbeta.obs
    # update A
    UDsq = UD(U,Dsq,n)
    M_obs = suvC(UDsq,V,irow,pcol)
    S@x = y@x - M_obs - xbeta.obs  
    A.partial = ((S%*%V) + UDsq)
    A = A.partial - svdH$u %*% (svdH$v %*% A.partial)
    if(lambda>0) A = t(t(A) * (Dsq/(Dsq+lambda))) 
    Asvd=  fast.svd(as.matrix(A))
    U = Asvd$u
    Dsq = Asvd$d
    V = V %*% (Asvd$v)
    HU = svdH$u %*% (svdH$v %*% U)
    #-----------------------------------------------------------------------------------
    if(trace.it | return_obj)  obj= (.5*sum(S@x^2)+lambda*sum(Dsq))/nz 
    if(return_obj) obj.l[iter] = obj
    #------------------------------------------------------------------------------
    ratio=  Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    #------------------------------------------------------------------------------
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio,"\n")
  
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if(lambda>0&final.svd){
  
    part1 = suvC(svdH$u, t(as.matrix(svdH$v %*% S)), irow, pcol)
    part2 = suvC(HU,UD(V,Dsq,m), irow, pcol)
    xbeta.obs = part1 + part2 + xbeta.obs
        
    UDsq = UD(U,Dsq,n)
    M_obs = suvC(UDsq,V,irow,pcol)
    S@x = y@x - M_obs - xbeta.obs  
    A.partial = ( (S%*%V) + UDsq )
    A = A.partial - svdH$u %*% (svdH$v %*% A.partial)
    Asvd=fast.svd(as.matrix(A))
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      UDsq=UD(U,Dsq,n)
      M_obs=suvC(UDsq,V,irow,pcol)
      S@x=y@x-M_obs - xbeta.obs
      obj = (.5*sum(S@x^2)+lambda*sum(Dsq))/nz
      if(return_obj) obj.l[iter+1] = obj
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
  }
  J=min(sum(Dsq>0)+1,J)
  J = min(J, length(Dsq))
  out=list(u=U[, seq(J), drop=FALSE], d=Dsq[seq(J)], v=V[,seq(J), drop=FALSE], lambda=lambda, J=J,
           n_iter=iter, xbeta.obs=xbeta.obs)
  if(return_obj) out$obj = obj.l 
  out
}
