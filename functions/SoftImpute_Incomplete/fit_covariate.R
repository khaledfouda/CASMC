require(corpcor)
simpute.als.fit_Incomplete <-
function (y, X, J = 2, thresh = 1e-05, lambda=0, 
          maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE) {
  ###This function expects an object of class "Incomplete" which inherits from "sparseMatrix", where the missing entries
  ###are replaced with zeros. If it was centered, then it carries the centering info with it
  
  #this.call=match.call()
  #a=names(attributes(x))
  #binames=c("biScale:row","biScale:column")
  #if(all(match(binames,a,FALSE))){
  #  biats=attributes(x)[binames]
  #} else biats=NULL
  
  timespent = rep(0,4)
  if(!inherits(y,"dgCMatrix")) y=as(y,"dgCMatrix")
  irow=y@i
  pcol=y@p
  n <- dim(y)
  m <- n[2]
  n <- n[1]
  nz=nnzero(y)
  #-------------------------------
  start_time <- Sys.time()
  # hat2 <- X %*% solve(t(X) %*% X) %*% t(X)
  # 
  #p = ncol(X)
  Q <- qr.Q(Matrix::qr(X)) #[,1:p]
  hat <- Q %*% t(Q)
  #svdX <- fast.svd(X)
  #hat = svdX$u %*% Diagonal(length(svdX$d)) %*% t(svdX$u)
  onesSparse <- ys
  onesSparse@x[] <- 1
  hat_sp <- as(as.matrix(hat), "Incomplete")
  hat_sp <- hat_sp * onesSparse
    
   
    
    timespent[1] = timespent[1] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
  
  # hat_sp[y == 0] = NA
  # hat_sp <- as(hat_sp, "Incomplete")
  #---------------------------------------------------
  warm=FALSE
  if(!is.null(warm.start)){
    clean.warm.start(warm.start)
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
      r = J
    }else{
      Dsq=c(D,rep(D[JD],J-JD))
      r = length(Dsq)
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
    # update yfill[ynas] = UDV + XBETA
  }else{
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
      # update yfill[ynas] = 0
    }
  
  r = length(Dsq)
  # initial beta estimates using yfill - yplus = yfill - xbeta
  ratio <- 1
  iter <- 0
  S=y
  
  #-------------
  
    start_time <- Sys.time()
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old=V
    Dsq.old=Dsq
    #---------------------------------
    ## U step # S is yplus
    if(iter>1|warm){
      # 1
      BD=UD(V,Dsq,m)
      # 2
      yfill=suvC(U,BD,irow,pcol)
      # 3
      S@x = y@x - yfill - lsq.sp
    }else BD=matrix(0,m,r) 
    # 4
    lsq.sp = hat_sp@x * S@x
    # 5
    #UtHat = t(U) %*% hat
    # 6
    #print(paste("Execution time 1 is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),5), "seconds"))
    tU = t(U)
    tBD = t(BD)
    B =    (tU - (tU %*% hat_sp)) %*% S         + 
      (tBD -  (((tU%*% hat) %*% U) %*% tBD) )
    #B=t(t(U)%*%S)+BD
    #print(paste("Execution time 2 is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),5), "seconds"))
    if(lambda>0)B= B * (Dsq/(Dsq+lambda)) #UD(B,Dsq/(Dsq+lambda),m)
    #print(paste("Execution time 3 is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),5), "seconds"))
    Bsvd=fast.svd((as.matrix(B)))
    V=(Bsvd$v)      #V=Bsvd$u
    Dsq=Bsvd$d
    U=U%*% t(Bsvd$u)
    r = length(Dsq)
    #print(paste("Execution time 4 is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),5), "seconds"))
    # update yhat = UDV then yfill[ynas] = yhat[ynas] + xbeta[ynas]
    # compute beta estimates and yplus (S) = yfill - xbeta
    #------------------------------------
    ## V step
    # 1
    AD = UD(U,Dsq,n)
    # 2
    yfill= suvC(AD,V,irow,pcol)
    # 3
    S@x = y@x - yfill - lsq.sp
  if(trace.it)  obj=(.5*sum(S@x^2)+lambda*sum(Dsq))/nz # update later
    # 4
    lsq.sp = hat_sp@x * S@x
    # 5
    SV = (S %*% V)
    #start_time <- Sys.time()
    A = (SV - hat_sp %*% SV) + (AD - hat %*% AD)
    # timespent[2] = timespent[2] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
    #A=S%*%V+AD
    Dsq_scaled <- sapply(Dsq, function(d) d / (d + lambda))
    
    if(lambda>0) A = t(t(A) * (Dsq/(Dsq+lambda))) #UD(A,Dsq/(Dsq+lambda),n)
    Asvd=  fast.svd(as.matrix(A))
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    r = length(Dsq)
    # update yhat = UDV then yfill[ynas] = yhat[ynas] + xbeta[ynas]
    # compute beta estimates and yplus (S) = yfill - xbeta
    #------------------------------------------------------------------------------
    ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  }
    timespent[3] = timespent[3] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
    start_time <- Sys.time()
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if(lambda>0&final.svd){
    AD=UD(U,Dsq,n)
    xfill=suvC(AD,V,irow,pcol)
    S@x=y@x-yfill #- lsq.sp
    A=S%*%V+AD
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      AD=UD(U,Dsq,n)
      yfill=suvC(AD,V,irow,pcol)
      S@x=y@x-yfill
      obj=(.5*sum(S@x^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
  }
  J=min(sum(Dsq>0)+1,J)
  out=list(u=U[, seq(J)], d=Dsq[seq(J)], v=V[,seq(J)], lambda=lambda, rank=J)
  #attributes(out)=c(attributes(out),list(lambda=lambda,call=this.call),biats)
  timespent[4] = timespent[4] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
  print(round(timespent, 3))
  out
}
