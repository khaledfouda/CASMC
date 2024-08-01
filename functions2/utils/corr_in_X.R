

utils$cor_in_X <- function(X, thresh = 0.5, ndig = 2) {
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

utils$remove_collinear_cols <-
  remove_collinear_cols <-
  function(matrix, thresh = 0.7) {
    cor_m <- cor(matrix)
    to_remove <-
      caret::findCorrelation(cor_m, thresh, verbose = TRUE)
    matrix <- matrix[, -to_remove, drop = FALSE]
    return(matrix)
  }