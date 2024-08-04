utils$inv <- function(X, is_square = nrow(X) == ncol(X)) {
  if (is_square) {
    # Try to apply solve() and catch any errors indicating singularity
    tryCatch({
      return(solve(X))
    }, error = function(e) {
      return(ginv(X))
    })
  } else {
    # Use ginv() for non-square matrices
    inv_X <-
      return(ginv(X))
  }
}


utils$prepare.M.warm.start <- function(warm.start, J, n, m, min_eigv) {
  D = warm.start$d
  JD = sum(D > min_eigv)
  if (JD >= J) {
    U = warm.start$u[, seq(J), drop = FALSE]
    V = warm.start$v[, seq(J), drop = FALSE]
    Dsq = D[seq(J)]
  } else{
    # upscale
    Ja = J - JD
    Dsq = c(D, rep(D[JD], Ja))
    U = warm.start$u
    Ua = matrix(rnorm(n * Ja), n, Ja)
    Ua = Ua - U %*% (t(U) %*% Ua)
    Ua = utils$svd_small_nc(Ua, trim=FALSE, p=Ja)$u
    U = cbind(U, Ua)
    V = cbind(warm.start$v, matrix(0, m, Ja))
  }
  list(U=U, Dsq=Dsq, V=V)
}






utils$Frob <-
  Frob <-
  function(Uold, Dsqold, Vold, U, Dsq, V) {
    denom = sum(Dsqold ^ 2)
    utu = Dsq * (t(U) %*% Uold)
    vtv = Dsqold * (t(Vold) %*% V)
    uvprod = sum(diag(utu %*% vtv))
    num = denom + sum(Dsq ^ 2) - 2 * uvprod
    num / max(denom, 1e-9)
  }


utils$reduced_hat_decomp <-
  reduced_hat_decomp <-
  function(X, tol = 1e-2, pct = NULL) {
    # returns X (in reduce rank), the SVD of the hat matrix, and the rank.
    # the pct is the percentage of singular values contribution in X.
    # if not specified, then no reduction is done, unless one of the
    # singular values is smaller than the tolerance value provided.
    
    X_svd = utils$fast.svd(X, tol)
    if (!is.null(pct)) {
      stopifnot(pct >= 0.5 && pct < 1)
      var_contrib = cumsum(X_svd$d / sum(X_svd$d))
      n_d <- sum(var_contrib < pct)
      X_svd$d <- X_svd$d[1:n_d]
      X_svd$u <- X_svd$u[, 1:n_d]
      X_svd$v <- X_svd$v[, 1:n_d]
    }
    X = X_svd$u %*% (X_svd$d * t(X_svd$v))
    rank = length(X_svd$d)
    print(paste("Rank of X changed from", min(dim(X)), "to", rank))
    Q <- qr.Q(Matrix::qr(X))
    H <- Q %*% t(Q)
    svdH <- tryCatch({
      irlba(H, nu = rank, nv = rank, tol = tol)
    },
    error = function(e) {
      message(paste("SvdH:", e))
      utils$svd_simple(H, rank)
    })
    svdH$u = svdH$u
    svdH$v = svdH$d * t(svdH$v)
    svdH$d <- H <- Q <- NULL
    return(list(X = X, svdH = svdH, rank = rank))
  }


utils$reduced_hat_decomp.H <-
  reduced_hat_decomp.H <-
  function(X) {
    # returns only the SVD of the hat matrix.
    qrX = Matrix::qr(X)
    rank = qrX$rank
    Q <- qr.Q(qrX)
    H <- Q %*% t(Q)
    svdH <- tryCatch({
      irlba(H, nu = rank, nv = rank, tol = 1e-5)
    },
    error = function(e) {
      message(paste("SvdH:", e))
      utils$svd_simple(H, rank)
    })
    list(u = svdH$u,
         v = svdH$d * t(svdH$v),
         rank = rank)
  }

utils$UD <- UD <- function(U, D, n = nrow(U)) {
  U * outer(rep(1, n), D, "*")
}


utils$GetXterms <- GetXterms <- function(X, lambda = 0) {
  svdX = utils$fast.svd(X)
  Ux = svdX$u
  Vx = svdX$d * t(svdX$v)
  X0 = ginv(t(Vx) %*% Vx + diag(lambda, ncol(X))) %*% t(Vx)
  X1 = X0 %*% t(Ux) # k by n
  X2 = X0 %*% Vx # k by k
  list(X1 = X1, X2 = X2)
}


utils$lambdaM.max <-
  function(Y,
           svdH = NULL,
           X = NULL,
           tol = 1e-1,
           max_it = 30) {
    # initial values for lambda_M when the Y is incomplete.
    if (!is.null(svdH)) {
      Y = Y - svdH$u %*% (svdH$v %*% Y)
    } else if (!is.null(X)) {
      Y = Y - X %*% utils$inv(t(X) %*% X, TRUE) %*% t(X) %*% Y
    }
    
    propack.svd(as.matrix(Y), 1, opts = list(kmax = max_it))$d[1]
    
  }

utils$lambda.beta.max <-
  function(Y,
           svdH = NULL,
           X = NULL,
           tol = 1e-1,
           max_it = 30) {
    # initial values for lambda_M when the Y is incomplete.
    if (!is.null(svdH)) {
      Xbeta =  svdH$u %*% (svdH$v %*% Y)
    } else if (!is.null(X)) {
      Xbeta =  X %*% utils$inv(t(X) %*% X, TRUE) %*% t(X) %*% Y
    }
    
    propack.svd(as.matrix(Xbeta), 1, opts = list(kmax = max_it))$d[1]
    
  }




suvC <- utils$suvC <-
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
