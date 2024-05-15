naive_fit <- function(y, X) {
 X_r <- reduced_hat_decomp(X)
 Y_naive = as.matrix(y)
 Y_naive = naive_MC(Y_naive)
 naive_fit <-  X_r$svdH$u %*% (X_r$svdH$v  %*% Y_naive)
 naive_fit <- Y_naive - naive_fit
 #----------------------
 # initialization for beta = X^-1 Y
 # comment for later: shouldn't be X^-1 H Y??
 beta = as.matrix(ginv(X) %*% Y_naive)
 return(list(
  estimates = Y_naive,
  M = naive_fit,
  beta = beta
 ))
}
