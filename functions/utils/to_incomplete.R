to_incomplete <- function(mat) {
 mat[mat == 0] = NA
 mat <- as(mat, "Incomplete")
}
unsvd <- function(l, make_sparse=FALSE){
 res <- l$u %*% (l$d * t(l$v))
 if(make_sparse){
  row_contrib = rowSums(l$u^2)
  sorted_ind = order(row_contrib, decreasing = TRUE)
  
  rows_to_remove = (length(l$d)+1):nrow(l$u)
  #rows_to_remove = which((cumsum(row_contrib[sorted_ind]) / sum(row_contrib)) > thresh)
  
  rows_to_remove = sorted_ind[rows_to_remove]
  res[rows_to_remove,] <- 0
  }
 return(res)
}
svd_trunc_simple <- function(mat, J){
  mat_svd = svd(mat)
  mat_svd$d <- mat_svd$d[1:J, drop=FALSE]
  mat_svd$u <- mat_svd$u[, 1:J, drop=FALSE]
  mat_svd$v <- mat_svd$v[, 1:J, drop=FALSE]
  return(mat_svd)
}
