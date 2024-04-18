reduced_hat_decomp <- function(X, tol = 1e-2) {
   X_svd = fast.svd(X, tol)
   X = X_svd$u %*% (X_svd$d * t(X_svd$v))
   rank = length(X_svd$d)
   Q <- qr.Q(Matrix::qr(X))
   H <- Q %*% t(Q)
   svdH <- fast.svd(H, tol)
   svdH$u = svdH$u
   svdH$v = svdH$d * t(svdH$v)
   svdH$d <- H <- Q <- NULL
   return(list(X = X, svdH = svdH, rank = rank))
}
UD = function(U, D, n = nrow(U)) {
   U * outer(rep(1, n), D, "*")
}
GetXterms <- function(X) {
   svdX = fast.svd(X)
   Ux = svdX$u
   Vx = svdX$d * t(svdX$v)
   X0 = ginv(t(Vx) %*% Vx) %*% t(Vx)
   X1 = X0 %*% t(Ux)
   X2 = X0 %*% Vx
   list(X1 = X1, X2 = X2)
}
lambda0.cov_splr <- function(Y,
                             svdH,
                             tol = 1e-1,
                             max_it = 30) {
   yplus = Y - svdH$u %*% (svdH$v %*% Y)
   if (inherits(Y, "dgCMatrix") & (dim(Y)[1] == dim(Y)[2]) ) {
      eigs_sym(yplus, k = 1, which = "LM")$values
   } else
      propack.svd(as.matrix(yplus), 1, opts = list(kmax = max_it))$d[1]
}




suvC <-
   function(u, v, irow, pcol) {
      dd = dim(u)
      nnrow = as.integer(dd[1])
      nncol = as.integer(nrow(v))
      nrank = dd[2]
      storage.mode(u) = "double"
      storage.mode(v) = "double"
      storage.mode(irow) = "integer"
      storage.mode(pcol) = "integer"
      nomega = as.integer(length(irow))
      .Fortran(
         "suvC",
         nnrow,
         nncol,
         nrank,
         u,
         v,
         irow,
         pcol,
         nomega,
         r = double(nomega),
         PACKAGE = "softImpute"
      )$r
   }



# sparse_prod <-
#   function(n,m,r,H,sp,si,sx){
#     storage.mode(H)="double"
#     storage.mode(sx)="double"
#     storage.mode(si)="integer"
#     storage.mode(sp)="integer"
#     storage.mode(n)="integer"
#     storage.mode(m)="integer"
#     storage.mode(r)="integer"
#     result = double(r)
#
#     .Fortran("sparse_prod",
#              n,m,r,H,sp,si,sx
#     )$result
#   }
#
#
# sparse_prod <- function(n, m, r, H, sp, si, sx) {
#   # Ensure correct data types
#   storage.mode(H) = "double"
#   storage.mode(sx) = "double"
#   storage.mode(si) = "integer"
#   storage.mode(sp) = "integer"
#   storage.mode(n) = "integer"
#   storage.mode(m) = "integer"
#   storage.mode(r) = "integer"
#
#   # Preallocate the result vector
#   result = double(r)
#
#   # Call the Fortran subroutine
#   .Fortran("sparse_prod",
#            as.integer(n), as.integer(m), as.integer(r),
#            as.double(H), as.integer(sp), as.integer(si),
#            as.double(sx), result = as.double(result)
#   )$result
# }
