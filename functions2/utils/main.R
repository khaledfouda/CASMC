utils <- new.env()




library(irlba)
library(RSpectra)
library(Matrix)
library(microbenchmark)



mat_thin <- matrix(rnorm(800*10), ncol = 10) 
mat_dense <- matrix(rnorm(800*900), ncol = 900)


mat <- m <- mat_dense
nv <- 10
nr <- 800
nc = 900

microbenchmark(
 base = utils$svd_simple(mat, nv),
 fast = fast.svd(mat),
 propack = propack.svd(mat, nv),
 #irl = irlba(mat, nv),
 ss = svds(mat, nv),
 new = utils$svdopt(mat, nv, nr, nc),
 times = 100
)

n = dim(m)[1]
p = dim(m)[2]


EDGE.RATIO = 2 # use standard SVD if matrix almost square
if (n > EDGE.RATIO * p)
{
 return(psmall.svd(m, tol, trim, p))
}
else if (EDGE.RATIO * n < p