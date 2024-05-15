to_incomplete <- function(mat) {
 mat[mat == 0] = NA
 mat <- as(mat, "Incomplete")
}
unsvd <- function(l){
 l$u %*% (l$d * t(l$v))
}
