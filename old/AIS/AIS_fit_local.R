require(Matrix)

AIS_fit <- function(Y, lambda, maxR = NA, maxIter=100, tol=1e-6,
                    thresh=1e-3,decay=0.8, trace.it=FALSE, exact=FALSE){
   
   yfill <- Y # D is Y
   ymiss <- Y == 0
   m <- dim(Y)[1]
   n <- dim(Y)[2]
   if(is.na(maxR)) maxR = min(dim(Y))
   
   #lambdaMax <- svd(Y, nu = 1, nv = 0)$d[1] # replace
   lambdaMax <- topksvd(Y, 1, 5) 
   
   # sparse matrix
   R <- matrix(rnorm(n), n, 1) 
   D_sparse <- as(Y, "sparseMatrix") # this is Z
   rows <- row(D_sparse)
   cols <- col(D_sparse)
   data <- D_sparse@x
   
   
   #Q <- powerMethod(Y, R, maxIter=5, tol=tol)
   U0 <- U1 <- powerMethod(Y, R, maxIter=5, tol=tol)
   V0 <- V1 <- svd(t(U0) %*% Y)$v # nu=0? why
   
   a0 <- a1 <- 1
   
   # M <- M.old <- matrix(0, m, n)
   # iter <- 0
   # ratio <- 1
   
   # Objectives and metrics
   obj <- numeric(maxIter)
   RMSE <- numeric(maxIter)
   Time <- numeric(maxIter)
   RankIn <- numeric(maxIter)
   RankOut <- numeric(maxIter)
   t <- Sys.time()
   # i =1
   for(i in 1:maxIter){
   #while((ratio > thresh)&(iter < maxit)){
     
     
      #svtY_old = svtY
      #iter <- iter + 1
      
      
      
      lambdai <- abs(lambdaMax - lambda) * (decay ^ i) + lambda
      bi <- (a0 - 1) / a1
      
      #theta = (iter-1)/(iter+2)
      
      # Sparse term Z = U*V' + spa = M + theta(M - M.old)
      part0 <- partXY(t(U0), t(V0), row(D_sparse), col(D_sparse), length(D_sparse@x))
      part1 <- partXY(t(U1), t(V1), row(D_sparse), col(D_sparse), length(D_sparse@x))
      part0 <- D_sparse@x - (1 + bi) * t(part1) + bi * t(part0)
      D_sparse@x <- as.vector(part0)
      #setSval(D_sparse, part0)
      #D_sparse <- sparseMatrix(i = rows, j = cols, x = part0, dims = c(m, n))
      
      #Z = M + theta * (M - M.old)
      #Z = Y - (1 + bi)  
      
      #yfill[ymiss] = M[ymiss]
      #-----------------------------------------------------------------
      # Compute R using V1 and V0
      # R <- V0 - V1 %*% t(V1) %*% V0
      # R <- colSums(R^2) > 0
      # R <- cbind(V1, V0[,R])
      # R <- R[, 1:min(ncol(R), maxR)]
      R <- filterBase(V1, V0, 1e-6)
      R <- R[, 1:min(ncol(R), maxR), drop=FALSE]
      RankIn[i] <- ncol(R)
      #--------------------------------------
      pwTol = max(1e-6, lambdaMax*(0.95^i));
      
      if(exact == TRUE){
        # Exact Method
        next;
        Ui_S_Vi <- matcompSVD(U1, V1, U0, V0, D_sparse, bi, ncol(R))
        Ui <- Ui_S_Vi$U
        S <- Ui_S_Vi$S
        Vi <- Ui_S_Vi$V
        
        S <- diag(S)
        nnzS <- sum(S > lambdai)
        Ui <- Ui[, 1:nnzS]
        Vi <- Vi[, 1:nnzS]
        S <- S[1:nnzS]
        S <- S - lambdai
        S <- diag(S)
        
        Ui <- Ui %*% S
        pwIter <- Inf
        
      }else{
        # Approximate method
        list_Q_pwIter <- powerMethodAccMatComp(U1, V1, U0, V0, D_sparse, bi, R, 10, pwTol)
        
        Q <- list_Q_pwIter$Q
        pwIter <- list_Q_pwIter$maxIter
        
        hZ <- (1 + bi) * (t(Q) %*% (U1)) %*% t(V1) - bi *
          (t(Q) %*% (U0)) %*% t(V0) + t(Q) %*% (D_sparse)
        
        Ui_S_Vi <- SVT(hZ, lambdai)
        Ui <- Ui_S_Vi$U
        S <- Ui_S_Vi$S
        Vi <- Ui_S_Vi$V
        
        Ui <- Q %*% (Ui %*% S)
        
      }
      
      
      #---------------------------------------
      # svtY <- SVT_Approx(yfill, R, lambdai, maxit)
      # M.old <- M
      # M <- svtY$u %*% (svtY$d * t(svtY$v))
      #------------------------------------------------------------
      # Update V1 and V0
      # V0 <- V1
      # V1 = svtY$v
      U0 <- U1
      U1 <- Ui
      V0 <- V1
      V1 <- Vi
      #-------------------
      RankOut[i] <- length(S)
      
      ai <- (1 + sqrt(1 + 4 * a0^2)) / 2
      a0 <- a1
      a1 <- ai
      
      objVal <- partXY(t(Ui), t(Vi), row(D_sparse), col(D_sparse), length(D_sparse@x))
      objVal <- 0.5 * sum((D_sparse@x - objVal)^2)
      objVal <- objVal + lambda * sum(S)
      obj[i] <- objVal
      
      # find another way for the ratio as U and U.old won't have the same dimension
      # if(iter != 1)
      #    ratio=Frob(svtY_old$u,svtY_old$d,svtY_old$v,svtY$u,svtY$d,svtY$v)
      # svtY_old = svtY
      
      if (i > 1) {
        delta <- obj[i - 1] - objVal
        print(sprintf('iter: %g; obj: %3f (dif: %3f); rank %g; lambda: %1f; power(iter %g, rank %g, tol %2f)', 
                    i, objVal, delta, sum(S !=0), lambdai, pwIter, ncol(R), pwTol))
        
        # Adaptive restart
        if (delta < 0) {
          a0 <- 1
          a1 <- 1
        }
      } else {
        print(sprintf('iter: %g; obj: %g; rank %g; lambda: %1f; power(iter %g, rank %g, tol %2f)', 
                    i, objVal, sum(S != 0), lambdai, pwIter, ncol(R), pwTol))
      }
      # Checking convergence
      if (i > 1 && abs(delta) < tol) {
        break
      }
      
        
  }
  svd_U1 <- svd(U1)
  U0 <- svd_U1$u
  S <- svd_U1$d
  V <- svd_U1$v %*% t(V1)
  
  list(U0 = U0, S = S, V = V,
       output = list(obj = obj, Rank = RankOut,
                      RMSE = RMSE, RankIn = RankIn,
                      RankOut = RankOut, Time = Time))

   
}
