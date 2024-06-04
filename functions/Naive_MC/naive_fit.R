naive_fit <- function(y, X, return_xbeta = FALSE) {
  svdH <- reduced_hat_decomp.H(X)
  Y_naive = as.matrix(y)
  Y_naive = naive_MC(Y_naive)
  
  Xbeta <-  svdH$u %*% (svdH$v  %*% Y_naive)
  if (return_xbeta)
    return(Xbeta)
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
