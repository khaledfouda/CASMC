utils$generate_random_similarity_matrix <- function(dim, seed = NULL) {
 if(! is.null(seed)) set.seed(seed)
 matrix <- matrix(0, dim, dim)
 matrix[upper.tri(matrix)] <- sample(c(0, 1), sum(upper.tri(matrix)), replace = TRUE)
 matrix <- matrix + t(matrix) # make it symmetric
 diag(matrix) <- 1
 return(matrix)
}


utils$computeLaplacian <- function(S, normalized = TRUE) {
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


# 
# update_chol <- function(chol_M, D_vector) {
#  diag_update <- diag(sqrt(as.vector(D_vector)))  # Square root for Cholesky update
#  return(chol_M + diag_update)  
# }