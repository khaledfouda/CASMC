AIS_fit <- function(Y, lambda, maxR = NA, maxit=100, tol=1e-6, thresh=1e-3,decay=0.8, trace.it=FALSE){
   
   yfill <- Y
   ymiss <- Y == 0
   m <- dim(Y)[1]
   n <- dim(Y)[2]
   R <- matrix(rnorm(n), ncol = 1)
   if(is.na(maxR)) maxR = min(dim(Y))
   
   lambdaMax <- svd(Y, nu = 1, nv = 0)$d[1] # replace
   lambdaMax <- topksvd(Y, 1, 5) 
   
   
   Q <- powerMethod(Y, R, maxIter=5, tol=tol)
   U0 <- U1 <- powerMethod(Y, R, maxIter=5, tol=tol)
   
   
   V0 <- V1 <- svd(t(Q) %*% Y)$v # nu=0? why
   M <- M.old <- matrix(0, m, n)
   iter <- 0
   ratio <- 1
   
   
   a0 <- a1 <- 1
   
   while((ratio > thresh)&(iter < maxit)){
      #svtY_old = svtY
      iter <- iter + 1
      
      
      
      lambdai <- abs(lambdaMax - lambda) * (decay ^ iter) + lambda
      bi <- (a0 - 1) / a1
      
      theta = (iter-1)/(iter+2)
      
      
      Z = M + theta * (M - M.old)
      Z = Y - (1 + bi)  
      
      yfill[ymiss] = M[ymiss]
      #-----------------------------------------------------------------
      # Compute R using V1 and V0
      R <- V0 - V1 %*% t(V1) %*% V0
      R <- colSums(R^2) > 0
      R <- cbind(V1, V0[,R])
      R <- R[, 1:min(ncol(R), maxR)]
      #---------------------------------------
      svtY <- SVT_Approx(yfill, R, lambdai, maxit)
      M.old <- M
      M <- svtY$u %*% (svtY$d * t(svtY$v))
      #------------------------------------------------------------
      # Update V1 and V0
      V0 <- V1
      V1 = svtY$v
      #-------------------
      # find another way for the ratio as U and U.old won't have the same dimension
      if(iter != 1)
         ratio=Frob(svtY_old$u,svtY_old$d,svtY_old$v,svtY$u,svtY$d,svtY$v)
      svtY_old = svtY
      
   }
   
}
length(V0)
dim(svtY$v)
svtY$u 