naive_fit <- function(y, X) {
 X_r <- reduced_hat_decomp(X)
 Y_naive = as.matrix(y)
 Y_naive = naive_MC(Y_naive)
 
 Xbeta <-  X_r$svdH$u %*% (X_r$svdH$v  %*% Y_naive)
 M <- Y_naive - Xbeta
 #----------------------
 # initialization for beta = X^-1 Y
 # comment for later: shouldn't be X^-1 H Y??
 beta = as.matrix(ginv(X) %*% Xbeta)
 return(list(
  estimates = Y_naive,
  M = M,
  beta = beta
 ))
}
