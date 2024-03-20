biScaleMatrix <- function(mat) {
 if (!is.matrix(mat)) {
  stop("Input must be a matrix.")
 }
 
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
 
 return(
  list(
   scaledMat = mat,
   rowMeans = rowMeansOrig,
   colMeans = colMeansOrig,
   rowSDs = rowSDsOrig,
   colSDs = colSDsOrig
  )
 )
}


revertBiScaledMatrix <- function(scaledMat, params) {
 if (!is.matrix(scaledMat)) {
  stop("Input must be a matrix.")
 }
 
 # Unscaling
 mat <- sweep(scaledMat, 1, params$rowSDs, FUN = "*")
 mat <- sweep(mat, 2, params$colSDs, FUN = "*")
 
 # Un-centering
 mat <- sweep(mat, 1, params$rowMeans, FUN = "+")
 mat <- sweep(mat, 2, params$colMeans, FUN = "+")
 
 return(mat)
}
