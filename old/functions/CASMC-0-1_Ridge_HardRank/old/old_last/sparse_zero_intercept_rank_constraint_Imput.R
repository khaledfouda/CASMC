SZIRCI <- function(Y,
                   X,
                   k,
                   thresh = 1e-5,
                   maxit = 100,
                   trace.it = TRUE,
                   warm.start = NULL,
                   final.trim = TRUE) {
   #' This functions solves the model Y = X beta and returns beta (B below)
   # Input: X = Ux Vx, B = D V;  ||Y-XB^T||;  Y: Partially observed mxn;
   # X is nxk and B is mxk; X is given.
   if (!inherits(Y, "dgCMatrix"))
      Y = as(Y, "dgCMatrix")
   
   irow = Y@i
   pcol = Y@p
   #ynas = is.na(Y)
   n <- dim(Y)
   m <- n[2]
   n <- n[1]
   
   if (trace.it)
      nz = nnzero(Y)
   
   if (is.null(warm.start)) {
      svdX = fast.svd(X)
      Ux = svdX$u
      Vx = svdX$d * t(svdX$v)
      X0 = ginv(t(Vx) %*% Vx) %*% t(Vx)
      X1 = X0 %*% t(Ux)
      X2 = X0 %*% Vx
      B = t(ginv(X) %*% naive_MC(as.matrix(Y))) # B = (X^-1 Y)'
      Bsvd = fast.svd(B)
      
   } else{
      X1 = warm.start$X1
      X2 = warm.start$X2
      Bsvd = warm.start$Bsvd
   }
   
   ratio <- Inf
   iter <- 0
   S <- Y
   while ((ratio > thresh) & (iter < maxit)) {
      iter <- iter + 1
      U.old = Bsvd$u
      V.old = Bsvd$v
      D.old = Bsvd$d
      UD = t(Bsvd$d * t(Bsvd$u))
      S@x = Y@x - suvC(X %*% Bsvd$v, UD, irow, pcol)
      #--------------------------------
      B = t(X1 %*% S + X2 %*% Bsvd$v %*% t(UD))
      Bsvd = fast.svd(as.matrix(B))
      #-----------------------------------------------------------------
      ratio =  Frob(U.old, D.old, V.old, Bsvd$u, Bsvd$d, Bsvd$v)
      #------------------------------------------------------------------------------
      if (trace.it) {
         obj = (.5 * sum(S@x ^ 2)) / nz
         cat(iter, ":", "obj", format(round(obj, 5)), "ratio", ratio, "\n")
      }
      
      #------------------------------------------------------------------------------
   }
   if (iter == maxit &
       trace.it)
      warning(paste("Convergence not achieved by", maxit, "iterations"))
   
   if (final.trim) {
      J = min(sum(Bsvd$d > 0) + 1, k)
      J = min(J, length(Bsvd$d))
      Bsvd$u = Bsvd$u[, seq(J), drop = FALSE]
      Bsvd$v = Bsvd$v[, seq(J), drop = FALSE]
      Bsvd$d = Bsvd$d[seq(J)]
   }
   Bsvd$n_iter = iter
   #out = #list(Bsvd=Bsvd)
   Bsvd
}
