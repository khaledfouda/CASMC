reduced_hat_decomp <- function(X, tol = 1e-2) {
   X_svd = fast.svd(X, tol)
   X = X_svd$u %*% (X_svd$d * t(X_svd$v))
   rank = length(X_svd$d)
   Q <- qr.Q(Matrix::qr(X))
   H <- Q %*% t(Q)
   #svdH <- fast.svd(H, tol)
   svdH <- irlba(H, nu = rank, nv = rank, tol = tol)
   svdH$u = svdH$u
   svdH$v = svdH$d * t(svdH$v)
   svdH$d <- H <- Q <- NULL
   return(list(X = X, svdH = svdH, rank = rank))
}
UD = function(U, D, n = nrow(U)) {
   U * outer(rep(1, n), D, "*")
}
GetXterms <- function(X, lambda=0) {
   svdX = fast.svd(X)
   Ux = svdX$u
   Vx = svdX$d * t(svdX$v)
   X0 = ginv(t(Vx) %*% Vx + diag(lambda, ncol(X))) %*% t(Vx)
   X1 = X0 %*% t(Ux) # k by n
   X2 = X0 %*% Vx # k by k
   list(X1 = X1, X2 = X2)
}
lambda0.cov_splr <- function(Y,
                             svdH,
                             tol = 1e-1,
                             max_it = 30) {
   yplus = Y - svdH$u %*% (svdH$v %*% Y)
   propack.svd(as.matrix(yplus), 1, opts = list(kmax = max_it))$d[1]
   # if (inherits(Y, "dgCMatrix") & (dim(Y)[1] == dim(Y)[2]) ) {
   #    eigs(yplus, k = 1, which = "LM")$values
   # } else
} 

computeLaplacian <- function(S, normalized = TRUE) {
   d <- rowSums(S)
   if (normalized) {
      d_inv_sqrt <- 1 / sqrt(d)
      L <- diag(1, nrow(S)) - (d_inv_sqrt %*% t(d_inv_sqrt)) * S
   } else {
      L <- diag(d) - S
   }
   L[L==0] <- NA
   L <- as(L, "Incomplete")
   return(L)
}
update_chol <- function(chol_M, D_vector) {
   diag_update <- diag(sqrt(as.vector(D_vector)))  # Square root for Cholesky update
   return(chol_M + diag_update)  
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


generate_similarity_matrix <- function(n, seed = NULL) {
   if(! is.null(seed)) set.seed(seed)
   matrix <- matrix(0, n, n)
   matrix[upper.tri(matrix)] <- sample(c(0, 1), sum(upper.tri(matrix)), replace = TRUE)
   matrix <- matrix + t(matrix) # make it symmetric
   diag(matrix) <- 1
   return(matrix)
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
