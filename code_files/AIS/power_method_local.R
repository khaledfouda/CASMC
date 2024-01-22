powerMethod <- function(A, R, maxIter, tol = 1e-5) {
   Y <- A %*% R
   Q <- qr.Q(qr(Y, LAPACK = TRUE), complete=FALSE)
   err <- Inf
   for (i in 1:maxIter) {
      Y <- A %*% (t(A) %*% Q)
      iQ <- qr.Q(qr(Y, LAPACK = TRUE), complete = FALSE)
      
      # err is the largest singular value of  the difference
      err <- norm(iQ[, 1] - Q[, 1], type = "2")
      
      Q <- iQ
      
      if (err < tol) {
         break
      }
   }
   
   #Q
   list(Q = Q, maxIter = i)
}

powerMethodAccMatComp <- function(U1, V1, U0, V0, spa, bi, R, maxIter, tol) {
  Y <- U1 %*% ((t(V1) %*% R) * (1 + bi))
  Y <- Y - U0 %*% ((t(V0) %*% R) * bi)
  Y <- Y + spa %*% R
  
  Q <- qr.Q(qr(Y), complete = FALSE)
  err <- numeric(maxIter)
  
  for (i in 1:maxIter) {
    AtQ <- (1 + bi) * (t(Q) %*% U1) %*% t(V1)
    AtQ <- AtQ - bi * (t(Q) %*% U0) %*% t(V0)
    AtQ <- t(AtQ)
    
    Y <- U1 %*% ((t(V1) %*% AtQ) * (1 + bi))
    Y <- Y - U0 %*% ((t(V0) %*% AtQ) * bi)
    Y <- Y + spa %*% AtQ
    
    iQ <- qr.Q(qr(Y), complete = FALSE)
    
    err[i] <- norm(iQ[, 1] - Q[, 1], type = "2")
    Q <- iQ
    
    if (err[i] < tol) {
      break
    }
  }
  
  list(Q = Q, maxIter = i)
}
