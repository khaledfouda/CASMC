require(corpcor)
simpute.als.fit_Incomplete_2 <-
function (y, X, yvalid, J = 2, thresh = 1e-05, lambda=0, 
          maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE,
          patience=3) {

  
  if(!inherits(y,"dgCMatrix")) y=as(y,"dgCMatrix")
  irow=y@i
  pcol=y@p
  n <- dim(y)
  m <- n[2]
  n <- n[1]
  nz=nnzero(y)
  #-------------------------------
  Q <- qr.Q(Matrix::qr(X)) #[,1:p]
  H <- Q %*% t(Q)
  I_H <- Diagonal(n) - H
  H <- NULL
  
  #---------------------------------------------------
  warm=FALSE
  clean.warm.start(warm.start)
  if(!is.null(warm.start)){
    #must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))
      stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=sum(D>0)
    lsq.sp <- warm.start$lsq.sp
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
    # update yfill[ynas] = UDV + XBETA
  }else{
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=fast.svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VDsq and AB'=U Dsq V^T
      # update yfill[ynas] = 0
    }
  tU = t(U)
  r = length(Dsq)
  # initial beta estimates using yfill - yplus = yfill - xbeta
  ratio <- 1
  iter <- 0
  S=y
  #----------------------------------------
  counter = 0
  best_score = Inf
  #----------------------------------------
  
    #start_time <- Sys.time()
  while ((ratio > thresh)&(iter<maxit)&(counter < patience)) {
    iter <- iter + 1
    U.old=U
    V.old= t(V)
    Dsq.old=Dsq
    #---------------------------------
    ## U step # S is yplus
    if(iter>1|warm){
      # 1
      VDsq=UD(V,Dsq,m)
      # 2
      M_obs = suvC(U,VDsq,irow,pcol)
      # 3
      S@x = y@x - M_obs
      
    }else VDsq=matrix(0,m,r) 
    # 6
    HU = I_H %*% U
    B = t(S) %*% HU + (VDsq %*% (tU %*% HU))
    
    if(lambda>0) B = t(t(B) * (Dsq/(Dsq+lambda))) 
    
    Bsvd=fast.svd((as.matrix(B)))
    V=(Bsvd$u)      #V=BsVDsq$u
    Dsq=Bsvd$d
    U=U%*% (Bsvd$v)
    r = length(Dsq)
    #------------------------------------
    ## V step
    # 1
    UDsq = UD(U,Dsq,n)
    # 2
    M_obs = suvC(UDsq,V,irow,pcol)
    # 3
    S@x = y@x - M_obs
  
    if(trace.it)  obj=(.5*sum(S@x^2)+lambda*sum(Dsq))/nz # update later
    # 4
    A = I_H %*% ( (S%*%V) + UDsq )
    
    if(lambda>0) A = t(t(A) * (Dsq/(Dsq+lambda))) 
    #-----------------------------------------------------------------------------------
    valid_preds = yvalid@x - suvC(UDsq, V, yvalid@i, yvalid@p)
    valid_error = sqrt(mean(valid_preds^2))
    if(valid_error < best_score){
      counter = 0
      best_score = valid_error
      best_iter = iter
    }else
      counter = counter + 1
    #-----------------------------------------------------------------------------------
    Asvd=  fast.svd(as.matrix(A))
    U= (Asvd$u)
    Dsq=Asvd$d
    V=V %*% (Asvd$v)
    r = length(Dsq)
    tU = t(U)
    #------------------------------------------------------------------------------
    ratio= 2#Frob2(U.old,Dsq.old,V.old,tU,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, 
                    "training RMSE",sqrt(mean((S@x)^2)),
                    "valid RMSE", valid_error, "\n")
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if(lambda>0&final.svd){
    UDsq=UD(U,Dsq,n)
    M_obs=suvC(UDsq,V,irow,pcol)
    S@x=y@x-M_obs #- lsq.sp
    A= S%*%V+UDsq
    Asvd=fast.svd(as.matrix(A))
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      UDsq=UD(U,Dsq,n)
      M_obs=suvC(UDsq,V,irow,pcol)
      S@x=y@x-M_obs
      obj=(.5*sum(S@x^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
  }
  J=min(sum(Dsq>0)+1,J)
  J = min(J, length(Dsq))
  out=list(u=U[, seq(J)], d=Dsq[seq(J)], v=V[,seq(J)], lambda=lambda, rank=J,
           lsq.sp = lsq.sp, best_iter=best_iter, best_score=best_score, last_score=valid_error)
  
  out
}
