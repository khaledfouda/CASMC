simpute.als.splr.fit.nocov.fixedJ <- function(X, J, thresh=1e-5, maxit=100, trace.it=TRUE,
                                              warm.start=NULL, final.trim=TRUE, return_obj=FALSE){
   if(!inherits(X, "dgCMatrix")) X = as(X, "dgCMatrix")
   
   irow = X@i
   pcol = X@p
   n <- dim(X)
   m <- n[2]
   n <- n[1]
   if(trace.it | return_obj) nz = nnzero(X)
   
   if(is.null(warm.start)){
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=fast.svd(U)$u
      Dsq=rep(1,J)
   }else{
      U = warm.start$u
      V = warm.start$v
      Dsq = warm.start$d
   }
   
   ratio <- Inf
   iter <- 0
   S = X
   if(return_obj) obj.l <- rep(NA, maxit)
   
   while((ratio > thresh)&(iter < maxit)){
      
      iter <- iter + 1
      U.old = U
      V.old = V
      Dsq.old = Dsq
      #--------------------------
      VDsq <- UD(V, Dsq, m)
      S@x = X@x - suvC(U, VDsq, irow, pcol)
      BD = as.matrix(t(S)%*%U + VDsq)
      #--------------------
      Bsvd=fast.svd(BD)
      V = Bsvd$u      
      Dsq = Bsvd$d
      U = U%*% (Bsvd$v)
      #----
      UDsq = UD(U, Dsq, n)
      S@x = X@x - suvC(UDsq, V, irow, pcol)
      AD = as.matrix(S%*%V + UDsq)
      #----
      Asvd=  fast.svd(AD)
      U = Asvd$u
      Dsq = Asvd$d
      V = V %*% (Asvd$v)
      #-----------------------------------------------------------------------------------
      ratio=  Frob(U.old,Dsq.old,V.old,U,Dsq,V)
      #------------------------------------------------------------------------------
      if(trace.it | return_obj)  obj= (.5*sum(S@x^2))/nz 
      if(return_obj) obj.l[iter] = obj
      if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio,"\n")
      #------------------------------------------------------------------------------
   }
   if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
   
   if(final.trim){
      J=min(sum(Dsq>0)+1,J)
      J = min(J, length(Dsq))
      U = U[,seq(J), drop=FALSE]
      V = V[,seq(J), drop=FALSE]
      Dsq = Dsq[seq(J)]
   }
   
   out = list(u=U, v=V, d=Dsq, J=J, iter=iter )
   if(return_obj) out$obj = obj.l
   out
}
