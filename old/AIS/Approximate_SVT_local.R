SVT_Approx <- function(A, R, lambda, maxIter=100, tol=1e-6) {
   
   
   # the following part computes Algorithm 1: Compute Q
   
   Q <- powerMethod(A, R, maxIter, tol)
   # 
   # Q <- qr(A %*% R, LAPACK = TRUE)$Q
   # 
   # err <- numeric(maxIter)
   # for (i in seq_len(maxIter)) {
   #    iQ <- qr(A %*% t(A %*% Q), LAPACK = TRUE)$Q
   #    
   #    err[i] <- norm(iQ[,1] - Q[,1], type = "2")
   #    Q <- iQ
   #    
   #    if(err[i] < tol) {
   #       break
   #    }
   # }
   #------------
   # the next part computes Algorithm 2:  Stage 2
   svdQA = svd(t(Q)%*%A)
   j = sum( (svdQA$d-lambda) > 0 )
   return(list(u = (Q%*%svdQA$u)[,1:j], d=svdQA$d[1:j], v=svdQA$v[,1:j]))
}

#A = yfill;lambda = lambdai; maxIter=100; tol=1e-6

topksvd <- function(X, k, maxIter = 5) {
  N <- ncol(X)
  R <- matrix(rnorm(N * k), N, k)
  
  Q <- powerMethod(X, R, maxIter)$Q
  
  X_Q <- t(Q) %*% X
  s <- svd(X_Q)$d
  return(s)
}
filterBase <- function(V1, V0, tol) {
  R <- V0 - V1 %*% (t(V1) %*% V0)
  R <- colSums(R^2)
  R <- R > tol
  R <- cbind(V1, V0[, R])
  return(R)
}

SVT <- function(Z, lambda, rnk = NULL) {
  if (!is.null(rnk)) {
    svd_Z <- svd(Z, nu = rnk, nv = rnk)
  } else {
    svd_Z <- svd(Z)
  }
  
  s <- svd_Z$d
  s <- pmax(s - lambda, 0)
  svs <- sum(s > 0)
  
  U <- svd_Z$u[, 1:svs, drop=FALSE]
  V <- svd_Z$v[, 1:svs, drop=FALSE]
  S <- s[1:svs] #diag(s[1:svs], svs)
  
  list(U = U, S = S, V = V)
}
partXY <- function(Xt, Y, I, J, L) {
  r <- nrow(Xt)
  m <- ncol(Xt)
  n <- ncol(Y)
  
  if (r != nrow(Y)) {
    stop("Rows of Xt must be equal to rows of Y")
  }
  if (r > m || r > n) {
    stop("Rank must be r <= min(m, n)")
  }
  
  v <- numeric(L)
  #V <- matrix(nrow = L, ncol = 1)
  
  for (p in 1:L) {
    ir <- (I[p] - 1) * r
    jr <- (J[p] - 1) * r
    #v[p] <- sum(Xt[ir + (1:r)] * Y[jr + (1:r)])
    v[p] <- sum(Xt[(ir+1):(ir+r)] * Y[(jr+1):(jr+r)])
    
  }
  
  v
}
