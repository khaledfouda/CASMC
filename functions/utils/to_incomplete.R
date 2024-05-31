to_incomplete <- function(mat) {
  mat[mat == 0] = NA
  mat <- as(mat, "Incomplete")
}
unsvd <- function(l, make_sparse = FALSE) {
  res <- l$u %*% (l$d * t(l$v))
  if (make_sparse) {
    row_contrib = rowSums(l$u ^ 2)
    sorted_ind = order(row_contrib, decreasing = TRUE)
    
    rows_to_remove = (length(l$d) + 1):nrow(l$u)
    #rows_to_remove = which((cumsum(row_contrib[sorted_ind]) / sum(row_contrib)) > thresh)
    
    rows_to_remove = sorted_ind[rows_to_remove]
    res[rows_to_remove,] <- 0
  }
  return(res)
}
svd_trunc_simple <- function(mat, J) {
  mat_svd = svd(mat)
  mat_svd$d <- mat_svd$d[1:J, drop = FALSE]
  mat_svd$u <- mat_svd$u[, 1:J, drop = FALSE]
  mat_svd$v <- mat_svd$v[, 1:J, drop = FALSE]
  return(mat_svd)
}



cor_in_cov <- function(X, thresh = 0.5, ndig = 2) {
  if (is.null(colnames(X)))
    colnames(X) <- paste0("Var", 1:ncol(X))
  cor_mat <- cor(X)
  cor_mat[lower.tri(cor_mat, diag = TRUE)] <- NA
  pairs <- which(abs(cor_mat) >= thresh, arr.ind = TRUE)
  data.frame(
    colname1 = colnames(X)[pairs[, 1]],
    colname2 = colnames(X)[pairs[, 2]],
    col_order1 = pairs[, 1],
    col_order2 = pairs[, 2],
    correlation = cor_mat[pairs] |> round(ndig)
  )
}


scalers <- function(X,
                    type = "median",
                    skip_bin = TRUE)
{
  is.mat = FALSE
  if (!is.null(ncol(X)) && ncol(X) > 1)
    is.mat = TRUE
  if (type == "median") {
    fun = function(x)
      (x - median(x)) / IQR(x)
  } else if (type == "minmax") {
    fun = function(x)
      (x - min(x)) / (max(x) - min(x))
  } else if (type == "standard") {
    fun = function(x)
      (x - mean(x)) / sd(x)
  } else if (type == "log") {
    fun = function(x)
      log1p(x + ifelse(min(x) < 0, abs(min(x)), 0))
  } else
    stop("type should be one of: median, minmax, standard or log ")
  
  if (is.mat == FALSE)
    if (length(unique(X)) <= 2 && skip_bin) {
      return(X)
    } else
      return(fun(X))
  
  fun2 = function(x, f, skip_bin) {
    if (length(unique(x)) <= 2 && skip_bin) {
      return(x)
    } else
      return(fun(x))
  }
  
  return(apply(X, 2, fun2, f = fun, skip_bin = skip_bin))
}
