lambda0.cov_splr <- function(Y,
                             svdH,
                             tol = 1e-1,
                             max_it = 30) {
 yplus = Y - svdH$u %*% (svdH$v %*% Y)
 if (inherits(Y, "dgCMatrix")) {
  eigs_sym(yplus, k = 1, which = "LM")$values
 } else
  propack.svd(yplus, 1, opts = list(kmax = max_it))$d[1]
}
