powerMethod <- function(A, R, maxIter, tol = 1e-5) {
   Y <- A %*% R
   Q <- qr.Q(qr(Y, LAPACK = TRUE))
   err <- Inf
   for (i in 1:maxIter) {
      Y <- A %*% (t(A) %*% Q)
      iQ <- qr.Q(qr(Y, LAPACK = TRUE))
      
      # err is the largest singular value of  the difference
      err <- norm(iQ[, 1] - Q[, 1], type = "2")
      
      Q <- iQ
      
      if (err < tol) {
         break
      }
   }
   
   #list(Q = Q, maxIter = i)
   Q
}
