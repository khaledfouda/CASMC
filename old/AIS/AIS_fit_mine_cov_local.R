require(Matrix)

AIS_cov_fit <- function(Y, X, lambda, maxR = NA, maxIter=100, tol=1e-6,
                    thresh=1e-3,decay=0.8, trace.it=FALSE, exact=FALSE){
   
   #Y <- Y # D is Y
   ymiss <- Y == 0
   m <- dim(Y)[1]
   n <- dim(Y)[2]
   if(is.na(maxR)) maxR = min(dim(Y))
   #-----------------------------------
   Y_x <- Y
   beta_partial = solve(t(X) %*% X) %*% t(X)
   beta.estim <- beta_partial %*% Y_x
   Xbeta <- X %*% beta.estim
   Y <- Y_x - Xbeta 
   
   #-----------------------------------
   lambdaMax <- topksvd(Y, 1, 5) 
   R <- matrix(rnorm(n), n, 1) 
   U0 <- U1 <- powerMethod(Y, R, maxIter=5, tol=tol)$Q
   V0 <- V1 <- svd(t(U0) %*% Y)$v # nu=0? why
   a0 <- a1 <- 1
   X1 <- U1 %*% t(V1)
   X0 <- U0 %*% t(V0)
   S0 <- S <-  c(0)
   #lambdaMax <- svd(Y, nu = 1, nv = 0)$d[1] # replace
   
   # sparse matrix
   #D_sparse <- as(Y, "sparseMatrix") # this is Z
   #rows <- row(D_sparse)
   #cols <- col(D_sparse)
   #data <- D_sparse@x
   
   
   #Q <- powerMethod(Y, R, maxIter=5, tol=tol)
   
   
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
      
      
      # lambdai is lambdat and lambdamax is lambda hat. decay is v.
     
      lambdai <- abs(lambdaMax - lambda) * (decay ^ (i-1) ) + lambda
      # a0 and a1 are c
      theta <- (a0 - 1) / a1
      
      # X1 is Xt
      Zt = (1+theta) * X1 - theta * X0
      Y_x[ymiss] = Zt[ymiss] + Xbeta[ymiss]
      beta.estim = beta_partial %*% Y_x
      Xbeta = X %*% beta.estim
      Y <- Y_x - Xbeta
      
      
      
      
      #-----------------------------------------------------------------
      # Compute R using V1 and V0
      #R <- filterBase(V1, V0, 1e-6)
      R <- V0 - V1 %*% (t(V1) %*% V0)
      R <- colSums(R^2) > 1e-6
      R <- cbind(V1, V0[, R])
      RankIn[i] <- ncol(R)
      if(RankIn[i] > maxR){
        R <- R[, 1:maxR, drop=FALSE]
        RankIn[i] <- maxR
      }
      #--------------------------------------
      pwTol = max(1e-6, lambdaMax*(0.95^i));
      
      if(exact == TRUE){
        # Exact Method
        next;
      }else{
        # Approximate method
        
        list_Q_pwIter <- powerMethod(Y, R, 10, pwTol)
        Q <- list_Q_pwIter$Q
        pwIter <- list_Q_pwIter$maxIter
        
        #----
        #Ui_S_Vi <- SVT(t(Q) %*% Y, lambdai)
        svd_Z <- svd(t(Q) %*% Y)
      
        s <- svd_Z$d - lambdai
        #s <- pmax(s - lambda, 0)
        svs <- sum(s > 0)
        
        Ui <- svd_Z$u[, 1:svs, drop=FALSE]
        Vi <- svd_Z$v[, 1:svs, drop=FALSE]
        S0 <- S
        S <- c(s[1:svs])
        Ui <- (Q %*% Ui) * S
        
        #----
        #Ui <- Ui_S_Vi$U
        #S0 <- S
        #S <- Ui_S_Vi$S
        #Vi <- Ui_S_Vi$V
        # print(dim(Q))
        # print(dim(Ui))
        # print(dim(diag(S)))
        # print(Ui_S_Vi$S)
        # X1 = Ui %*% Vi = U1 V1 and X0 = U0 V0
      }
      
      #------------------------------------------------------------
      # Update V1 and V0
      U0 <- U1
      U1 <- Ui
      V0 <- V1
      V1 <- Vi
      X0 <- X1
      X1 <- U1 %*% t(V1)
      #X0 <- U0 %*% t(V0)
      #-------------------
      RankOut[i] <- length(S)
      
      ai <- (1 + sqrt(1 + 4 * a0^2)) / 2
      a0 <- a1
      a1 <- ai
      
      # you need to implement the objective!!!!
      objVal <- i#0.5 * sum((Y - X1)^2) + lambda * sum(abs(S))
      #objVal <- Frob(U0, S0, V0, U1, S, V1)
      obj[i] <- objVal
      
      # objVal <- partXY(t(Ui), t(Vi), row(D_sparse), col(D_sparse), length(D_sparse@x))
      # objVal <- 0.5 * sum((D_sparse@x - objVal)^2)
      # objVal <- objVal + lambda * sum(S)
      # obj[i] <- objVal
      
      # find another way for the ratio as U and U.old won't have the same dimension
      # if(iter != 1)
      #    ratio=Frob(svtY_old$u,svtY_old$d,svtY_old$v,svtY$u,svtY$d,svtY$v)
      # svtY_old = svtY
      
      if (i > 1) {
        delta <- obj[i - 1] - objVal
        if(trace.it)
          print(sprintf('iter: %g; obj: %3f (dif: %3f); rank %g; lambda: %1f; power(iter %g, rank %g, tol %2f)', 
                    i, objVal, delta, sum(S !=0), lambdai, pwIter, ncol(R), pwTol))
        
        # Adaptive restart
        if (delta < 0) {
          a0 <- 1
          a1 <- 1
        }
      } else {
        if(trace.it)
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
# 
# dim(V0)
# S
# S0
