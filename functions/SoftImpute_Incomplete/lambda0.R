lambda0.cov_splr <- function(Y, svdH, tol=1e-1){
   yplus = Y - svdH$u %*% (svdH$v %*% Y)
   fast.svd(yplus)$d[1]
}

library(RSpectra)
result <- eigs_sym(matrix, k = 1, which = "LM")

# Display the largest eigenvalue
result$values