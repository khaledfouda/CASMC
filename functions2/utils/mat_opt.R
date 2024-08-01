utils$to_incomplete <- function(mat) {
  stopifnot(is.matrix(mat))
  mat[mat == 0] = NA
  mat <- as(mat, "Incomplete")
}

utils$unsvd <- function(obj,
                        ub = FALSE,
                        make_sparse = FALSE) {
  if (ub) {
    # this to match the SVD decomposition of Beta in the nuclear norm case.
    obj$u <- obj$ub
    obj$v <- obj$vb
    obj$d <- obj$db ^ 2
  }
  mat <- obj$u %*% (l$d * t(obj$v))
  if (make_sparse) {
    # delete this part later.
    row_contrib = rowSums(obj$u ^ 2)
    sorted_ind = order(row_contrib, decreasing = TRUE)
    rows_to_remove = (length(obj$d) + 1):nrow(obj$u)
    rows_to_remove = sorted_ind[rows_to_remove]
    mat[rows_to_remove, ] <- 0
  }
  return(mat)
}


utils$svd_simple <- function(mat, rank) {
  mat_svd = svd(mat)
  mat_svd$d <- mat_svd$d[1:rank, drop = FALSE]
  mat_svd$u <- mat_svd$u[, 1:rank, drop = FALSE]
  mat_svd$v <- mat_svd$v[, 1:rank, drop = FALSE]
  return(mat_svd)
}
