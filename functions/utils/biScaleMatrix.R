biScaleMatrix <-
   function(mat,
            with_min_max = FALSE,
            min_max_only = FALSE) {
      if (!is.matrix(mat) && !inherits(mat, "dgCMatrix")) {
         stop("Input must be a matrix or Incomplete matrix.")
      }
      
      rowMeansOrig <-
         colMeansOrig <- rowSDsOrig <- colSDsOrig <- NULL
      if (!min_max_only) {
         # Handling for sparse matrices
         if (inherits(mat, "dgCMatrix")) {
            stop("No implemented yet.")
            
         } else {
            mat[mat == 0] <- NA
            
            # Centering
            rowMeansOrig <- rowMeans(mat, na.rm = TRUE)
            mat <- sweep(mat, 1, rowMeansOrig, FUN = "-")
            
            colMeansOrig <- colMeans(mat, na.rm = TRUE)
            mat <- sweep(mat, 2, colMeansOrig, FUN = "-")
            
            # Scaling
            rowSDsOrig <- apply(mat, 1, sd, na.rm = TRUE)
            rowSDsOrig[rowSDsOrig == 0] <- 1
            mat <- sweep(mat, 1, rowSDsOrig, FUN = "/")
            
            colSDsOrig <- apply(mat, 2, sd, na.rm = TRUE)
            colSDsOrig[colSDsOrig == 0] <- 1
            mat <- sweep(mat, 2, colSDsOrig, FUN = "/")
         }
         
      }
      minVal <- maxVal <- NULL
      if (with_min_max | min_max_only) {
         mat[mat == 0] <- NA
         minVal <- min(mat, na.rm = TRUE)
         maxVal <- max(mat, na.rm = TRUE)
         
         # Scale to (0,1)
         mat <- (mat - minVal) / (maxVal - minVal)
         
         # Ensure all values are strictly greater than 0 by adding a small constant
         mat <- mat + .Machine$double.eps
      }
      
      return(
         list(
            scaledMat = mat,
            rowMeans = rowMeansOrig,
            colMeans = colMeansOrig,
            rowSDs = rowSDsOrig,
            colSDs = colSDsOrig,
            minVal = minVal,
            maxVal = maxVal
         )
      )
   }





revertBiScaledMatrix <- function(scaledMat, params) {
   if (!is.matrix(scaledMat)) {
      stop("Input must be a matrix.")
   }
   if (!is.null(params$minVal)) {
      # Revert min-max scaling
      mat <-
         scaledMat - .Machine$double.eps  # Remove the small constant added earlier
      mat <- (mat * (params$maxVal - params$minVal)) + params$minVal
   }
   # Unscaling
   if (!is.null(params$rowSDs)) {
      mat <- sweep(scaledMat, 1, params$rowSDs, FUN = "*")
      mat <- sweep(mat, 2, params$colSDs, FUN = "*")
      
      # Un-centering
      mat <- sweep(mat, 1, params$rowMeans, FUN = "+")
      mat <- sweep(mat, 2, params$colMeans, FUN = "+")
   }
   
   return(mat)
}
