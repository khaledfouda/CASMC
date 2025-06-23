lambda0.cov <- function(Y, X, maxiter = 30) {
 Xbeta <- X %*% ginv(t(X) %*% X) %*% t(X) %*% Y
 yplus <- Y - Xbeta
 return(propack.svd(yplus, 1, opts = list(maxiter = maxiter))$d[1])
}