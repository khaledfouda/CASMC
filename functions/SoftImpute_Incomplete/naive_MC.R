naive_MC <- function(mat) {
   # Input assumes that empty cells are filled with NA so they don't affect the average
   # Calculate row and column means, excluding NAs
   row_means <- rowMeans(mat, na.rm = TRUE)
   col_means <- colMeans(mat, na.rm = TRUE)
   
   # Expand row and column means to match the matrix dimensions
   expanded_row_means <- matrix(row_means[rep(1:length(row_means), each = ncol(mat))], nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
   expanded_col_means <- matrix(col_means[rep(1:length(col_means), times = nrow(mat))], nrow = nrow(mat), ncol = ncol(mat))
   avg_means <- (expanded_row_means + expanded_col_means) / 2
   
   na_positions <- is.na(mat)
   mat[na_positions] <- avg_means[na_positions]
   
   return(mat)
}
