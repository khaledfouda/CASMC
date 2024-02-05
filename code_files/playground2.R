library(Rcpp)

Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector sparse_prod_cpp(NumericMatrix H, IntegerVector sp, IntegerVector si, NumericVector sx, int n, int m, int r) {
  NumericVector result(r);
  int index = 0;

  for(int j = 0; j < m; ++j) {
    int jstart = sp[j];
    int jend = sp[j + 1];

    if(jstart >= jend) continue;

    for(int ind = jstart; ind < jend; ++ind) {
      int si_sub = si[ind];
      double sx_sub = sx[ind];
      double sum = 0;

      for(int k = 0; k < n; ++k) {
        sum += sx_sub * H(k, si_sub);
      }

      result[index++] = sum;
    }
  }

  return result;
}
')


setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

n <- 500
m <- 500
gen.dat <- generate_simulation_data_ysf(1,n,m,5,10, missing_prob = 0.9,coll=F,seed=3023)
X <- gen.dat$X
H <- X %*% solve(t(X)%*% X) %*% t(X)
#H <- matrix(rnorm(n * n), n, n)
mask <- Matrix(rbinom(n * m, 1, 0.1), n, m)
S <- matrix(rnorm(n * m), n, m) * mask
S <- Matrix(S, sparse = TRUE)

sp <- S@p
si <- S@i
sx <- S@x
r <- length(sx)

# Call the C++ function
result_cpp <- sparse_prod_cpp(H, sp, si, sx, n, m, r)


result_cpp %>% length()
length(result)

all(round(result_cpp,5) == round(result0,5))
result0[1:5]
result_cpp[1:5]

system.time({
   for(i in 1:30) 
   result0 = (H %*% S)[S!=0]
})

system.time({
  for(i in 1:30) 
    result0 = (p1 %*% (t(p2) %*% S))[S!=0]
})


system.time({
   for(i in 1:30) 
      # result_cpp <- sparse_prod_cpp(H, sp, si, sx, n, m, r)
   result = suvC(p1, as.matrix(t(S)%*%p2), si, sp )
})


pca_result <- prcomp(H)


pca_result$sdev


pca_res <- prcomp(H, scale. = TRUE)
summ <- summary(pca_res)

k = sum(summ$importance > 1e-2)
H_reduced <- pca_result$x[, 1:k] 

dim(pca_result$x)


pca_result <- prcomp(H, scale. = F)
summ <- summary(pca_result)
selected_comps <- which(summ$importance[2, ] > 1e4)
selected_comps <- which(summ$importance[3,] <=.99)
H_reduced <- pca_result$x[, selected_comps] %*% t(pca_result$rotation[, selected_comps])
dim(H_reduced)
p1 = pca_result$x[, selected_comps]
p2 = (pca_result$rotation[, selected_comps])

dim(p2)

result = suvC(p1, as.matrix(t(S)%*%p2), si, sp )
result[1:10]
result0[1:10]
length(result0)



summ$importance %>% t()  %>% as.data.frame() %>% arrange(desc(`Proportion of Variance`)) %>% head(10)

svd(H)$d -> dd; dd[dd>1e-3]
#########################################










S = ys
S@x = fits$xbeta.obs
W = Y_train == 0


n <- nrow(X)
k <- ncol(X)
m <- ncol(S)

loss_function <- function(m_flat) {
  M_est <- matrix(m_flat, nrow = k, ncol = m)
  S_pred <- (X %*% M_est) * W
  loss <- sum((S_pred - S)^2)
  return(loss)
}

initial_guess <- runif(k * m)
optim_result <- optim(initial_guess, loss_function, method = "L-BFGS-B")

M_est <- matrix(optim_result$par, nrow = k, ncol = m)



preds = best_preds + M
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))








##################################################