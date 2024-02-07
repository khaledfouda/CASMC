reduced_hat_decomp <- function(X, tol=1e-2){
   
   X_svd = fast.svd(X,tol)
   X = X_svd$u %*% (X_svd$d * t(X_svd$v))
   
   Q <- qr.Q(Matrix::qr(X)) 
   H <- Q %*% t(Q)
   svdH <- fast.svd(H, tol)
   svdH$u = svdH$u
   svdH$v = svdH$d * t(svdH$v)
   svdH$d <- H <- Q <- NULL
   return(list(X=X, svdH=svdH))
}
