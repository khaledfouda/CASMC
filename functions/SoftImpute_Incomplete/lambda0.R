library(RSpectra)
lambda0.cov_splr <- function(Y, svdH, tol=1e-1){
   yplus = Y - svdH$u %*% (svdH$v %*% Y)
   eigs_sym(yplus, k = 1, which = "LM")$values
}
