simpute.als.splr.fit.beta <- function(Y,X, J, thresh=1e-5, maxit=100, trace.it=TRUE,
                                              warm.start=NULL,
                                      final.trim=TRUE, return_obj=FALSE){
   # Input: X = Ux Vx, B = D V;  ||Y-XB^T||;  Y: Partially observed mxn; 
   # X is nxk and B is mxk; X is given.
   
   ynas = is.na(Y)
   n <- dim(Y)
   m <- n[2]
   n <- n[1]
   k <- J
   
   if(trace.it) nz = length(ynas)#nnzero(X)
   
   if(is.null(warm.start)){
      Y = naive_MC(Y)
      svdX = propack.svd(X,9)#fast.svd(X)
      Ux = svdX$u
      Vx = svdX$d * t(svdX$v)
      B = (ginv(t(Vx) %*% Vx) %*% t(Vx)) %*% t(Ux) %*% Y
      #B = t(ginv(X) %*% Y)
      #svdB = fast.svd(B)
      print(dim(B))
      svdB = propack.svd(t(B), 10)
      Vx = Vx %*% svdB$v
      D = svdB$d
      V = svdB$u
      
   }else{
      Ux = warm.start$Ux
      # complete later
   }
   
   ratio <- Inf
   iter <- 0

   while((ratio > thresh)&(iter < maxit)){
         
      iter <- iter + 1
      U.old = Ux %*% Vx
      V.old = V
      D.old = D
      
      #print(dim(V.old))
      #--------------------------------
      B = (ginv(t(Vx) %*% Vx) %*% t(Vx)) %*% t(Ux) %*% Y
      #print(dim(B))
      svdB = propack.svd(t(B), 10) #fast.svd(t(B))
      Vx = Vx %*% svdB$v
      D = svdB$d
      V = svdB$u
      Yhat =   Ux %*% Vx %*% (D * t(V))
      Y[ynas] = Yhat[ynas]
      #print(dim(svdB$v))
      #-----------------------------------------------------------------
      ratio=  Frob(U.old,D.old,V.old,Ux %*% Vx,D,V)
      #------------------------------------------------------------------------------
      if(trace.it)  obj= 2#(.5*sum(S@x^2))/nz 
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
   
   out = list(ux=Ux, vx=Vx, d=D, v=V, Yhat=Y)
   out
}
