require(corpcor)
simpute.als.fit_splr <-
function (y, yvalid, X=NULL, H=NULL, J = 2, thresh = 1e-05, lambda=0, 
          maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE,
          patience=3) {

  start_time <- Sys.time() #<<<<<<<<<<<<<
  if(!inherits(y,"dgCMatrix")) y=as(y,"dgCMatrix")
  irow=y@i
  pcol=y@p
  n <- dim(y)
  m <- n[2]
  n <- n[1]
  nz=nnzero(y)
  yobs = y != 0
  #-------------------------------
  if(is.null(H)){
    stopifnot(!is.null(X))
    Q <- qr.Q(Matrix::qr(X)) #[,1:p]
    H <- Q %*% t(Q)
    Q <- NULL
    #H = X %*% solve(t(X) %*% X) %*% t(X)
  }
  I_H <- Diagonal(n) - H
  
  #H #<- NULL
  
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
  }
  # initial beta estimates using yfill - yplus = yfill - xbeta
  ratio <- 1
  iter <- 0
  S=y
  xbeta.obs =  (H %*% y)[yobs]
  #----------------------------------------
  counter = 0
  best_score = Inf
  best_iter = NA
  print(paste("Execution time 1 is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
  #----------------------------------------
  time2 = 0
  time3 = 0
  #&(counter < patience)
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old= V
    Dsq.old=Dsq
    #---------------------------------
    HU = H %*% U
    ## U step # S is yplus
    if(iter>1|warm){
      # 1
      VDsq=UD(V,Dsq,m)
      # 2
      M_obs = suvC(U,VDsq,irow,pcol)
      S@x = y@x - M_obs - xbeta.obs  
  
    }else {
      VDsq=matrix(0,m,J) 
      S@x = y@x - xbeta.obs
    }
    # 6
    part1 <- numeric(length(yobs))
    start_time <- Sys.time() #<<<<<<<<<<<<
    
    for (j in 1:ncol(S)) {
      col_product <- H %*% S[, j]
      part1[which(S[, j] != 0)] <- col_product[yobs[, j]]
    }
    # part1 = suvC(H) (H %*% (S))[yobs] 
    time2 = time2 + round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2)
    start_time <- Sys.time() #<<<<
    part2 = suvC(HU,VDsq, irow, pcol)
    time3 = time3 + round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2)
    xbeta.obs = part1+part2 +xbeta.obs
    # xbeta.obs = (H %*% (S))[yobs] + suvC(HU,VDsq, irow, pcol) + xbeta.obs
    
    IHU =  U - HU
    
    B = as.matrix(t(S) %*% IHU) + (VDsq %*% (t(U) %*% IHU))  
    
    if(lambda>0) B = t(t(B) *(Dsq/(Dsq+lambda))) 
    
    Bsvd=fast.svd(B)
    V = Bsvd$u      
    Dsq = Bsvd$d
    U = U%*% (Bsvd$v)
    #------------------------------------
    ## V step
    
    
    UDsq = UD(U,Dsq,n)
    M_obs = suvC(UDsq,V,irow,pcol)
    S@x = y@x - M_obs - xbeta.obs  
    A = I_H %*% ( (S%*%V) + UDsq )
    if(lambda>0) A = t(t(A) * (Dsq/(Dsq+lambda))) 
    Asvd=  fast.svd(as.matrix(A))
    
    
    #-----------------------------------------------------------------------------------
    if(trace.it)  obj=2#(.5*sum(S@x^2)+lambda*sum(Dsq))/nz # update later
    #valid_preds = yvalid@x - suvC(UDsq, V, yvalid@i, yvalid@p)
    valid_error = Inf#sqrt(mean(valid_preds^2))
    if(valid_error < best_score){
      counter = 0
      best_score = valid_error
      best_iter = iter
    }else
      counter = counter + 1
    #-----------------------------------------------------------------------------------
    U = Asvd$u
    Dsq = Asvd$d
    V = V %*% (Asvd$v)
    
    #------------------------------------------------------------------------------
    ratio= 2# Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, 
                    #"training RMSE",sqrt(mean((S@x)^2)),
                    "valid RMSE", valid_error, "\n")
  
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if(lambda>0&final.svd){
    
    UDsq = UD(U,Dsq,n)
    M_obs = suvC(UDsq,V,irow,pcol)
    S@x = y@x - M_obs - xbeta.obs  
    A = I_H %*% ( (S%*%V) + UDsq )
    
    # UDsq=UD(U,Dsq,n)
    # M_obs=suvC(UDsq,V,irow,pcol)
    # S@x=y@x-M_obs #- HM_obs_sum_A
    # A= S%*%V+UDsq
    
    
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
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
  }
  
  print(paste("Execution time 2 is",time2, "seconds"))
  print(paste("Execution time 3 is",time3, "seconds"))
  
  J=min(sum(Dsq>0)+1,J)
  J = min(J, length(Dsq))
  out=list(u=U[, seq(J)], d=Dsq[seq(J)], v=V[,seq(J)], lambda=lambda, rank=J,
           best_iter=best_iter, best_score=best_score, last_score=valid_error,
           xbeta.obs=xbeta.obs)
  
  out
}
