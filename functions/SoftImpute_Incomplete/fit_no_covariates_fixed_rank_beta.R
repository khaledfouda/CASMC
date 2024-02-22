simpute.als.splr.fit.beta <- function(K,X, Xinv, J, thresh=1e-5, maxit=100, trace.it=TRUE,
                                              warm.start=NULL,
                                      final.trim=TRUE, return_obj=FALSE){
   knas = is.na(K)
   n <- dim(K)
   m <- n[2]
   n <- n[1]
   if(trace.it | return_obj) nz = length(knas)#nnzero(X)
   
   if(is.null(warm.start)){
      K = naive_MC(K)
      
   }else{
      K = warm.start
   }
   
   ratio <- Inf
   iter <- 0

   while((ratio > thresh)&(iter < maxit)){
      if(iter != 0){
         
      U.old = U
      V.old = V
      Dsq.old = Dsq
      }
      iter <- iter + 1
      #--------------------------
      ksvd = fast.svd(Xinv %*% K)
      U = ksvd$u
      V = ksvd$v
      Dsq = ksvd$d
      K[knas] <- ((X %*% U) %*% (Dsq * t(V) ))[knas]

      #-----------------------------------------------------------------------------------
      if(iter != 1){
         
         ratio=  Frob(U.old,Dsq.old,V.old,U,Dsq,V)
         #------------------------------------------------------------------------------
         if(trace.it)  obj= 2#(.5*sum(S@x^2))/nz 
         if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio,"\n")
      }
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
   
   out = list(u=U, v=V, d=Dsq, K=K, iter=iter )
   out
}
