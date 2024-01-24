simpute.als.fit_Incomplete <-
function (y, X, J = 2, thresh = 1e-05, lambda=0, 
          maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE) 
{
  ###This function expects an object of class "Incomplete" which inherits from "sparseMatrix", where the missing entries
  ###are replaced with zeros. If it was centered, then it carries the centering info with it
  
  #this.call=match.call()
  #a=names(attributes(x))
  #binames=c("biScale:row","biScale:column")
  #if(all(match(binames,a,FALSE))){
  #  biats=attributes(x)[binames]
  #} else biats=NULL
  
  if(!inherits(y,"dgCMatrix")) y=as(y,"dgCMatrix")
  irow=y@i
  pcol=y@p
  n <- dim(y)
  m <- n[2]
  n <- n[1]
  nz=nnzero(y)
  
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
    }
    else{
      Dsq=c(D,rep(D[JD],J-JD))
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
    # update yfill[ynas] = UDV + XBETA
  }
  else
    {
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
      # update yfill[ynas] = 0
    }
  
  # initial beta estimates using yfill - yplus = yfill - xbeta
  ratio <- 1
  iter <- 0
  yres=y
  
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old=V
    Dsq.old=Dsq
    #---------------------------------
    ## U step # yres is yplus
    if(iter>1|warm){
      BD=UD(V,Dsq,m)
      yfill=suvC(U,BD,irow,pcol)
      yres@x = y@x - yfill
    }else BD=0 
    B=t(t(U)%*%yres)+BD
    if(lambda>0)B=UD(B,Dsq/(Dsq+lambda),m)
    Bsvd=svd(B)
    V=Bsvd$u
    Dsq=Bsvd$d
    U=U%*%Bsvd$v
    # update yhat = UDV then yfill[ynas] = yhat[ynas] + xbeta[ynas]
    # compute beta estimates and yplus (yres) = yfill - xbeta
    #------------------------------------
    ## V step
    AD=UD(U,Dsq,n)
    yfill=suvC(AD,V,irow,pcol)
    yres@x = y@x - yfill
  if(trace.it)  obj=(.5*sum(yres@x^2)+lambda*sum(Dsq))/nz # update later
    A=yres%*%V+AD
    if(lambda>0)A= UD(A,Dsq/(Dsq+lambda),n)
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    # update yhat = UDV then yfill[ynas] = yhat[ynas] + xbeta[ynas]
    # compute beta estimates and yplus (yres) = yfill - xbeta
    #------------------------------------------------------------------------------
    ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if(lambda>0&final.svd){
    AD=UD(U,Dsq,n)
    xfill=suvC(AD,V,irow,pcol)
    yres@x=y@x-yfill
    A=yres%*%V+AD
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
    if(trace.it){
      AD=UD(U,Dsq,n)
      yfill=suvC(AD,V,irow,pcol)
      yres@x=y@x-yfill
      obj=(.5*sum(yres@x^2)+lambda*sum(Dsq))/nz
      cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }
  }
  J=min(sum(Dsq>0)+1,J)
  out=list(u=U[, seq(J)], d=Dsq[seq(J)], v=V[,seq(J)], lambda=lambda, rank=J)
  #attributes(out)=c(attributes(out),list(lambda=lambda,call=this.call),biats)
  out
}
