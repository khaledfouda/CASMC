

rm(list = ls())
library(softImpute)
library(purrr)
source("SMC_functions6.R")

set.seed(123)
n_rows <- 300
n_cols <- 300


X_cols <- 5
Z_cols <- 10

missing_prob <- 0.7
W_data <- matrix( rbinom(n_rows*n_cols, 1, (1 - missing_prob) ) , nrow = n_rows)

### sim1 data
#**************

sink("sim1_07maoxz.txt", append = TRUE)
print(c("mao_error_coop", "mao_error_x", "mao_error_z", "soft_error"))
sink()

for(iter_i in 1:50){
  
  print(iter_i)
  
  X_data <- matrix(rnorm(n_rows*X_cols), ncol = X_cols)
  Z_data <- matrix(rnorm(n_cols*Z_cols), ncol = Z_cols)
  E_data <- matrix(rnorm(n_rows*n_cols), ncol = n_cols)
  A <- matrix(rnorm(n_rows*10), nrow = n_rows)
  B <- matrix(rnorm(n_cols*10), ncol = n_cols)
  C <- matrix(rnorm(n_rows*10), nrow = n_rows)
  D <- matrix(rnorm(n_cols*10), ncol = n_cols)
  Px <- X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  Pz <- Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  Pxp <- diag(dim(X_data)[1]) - X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  Pzp <- diag(dim(Z_data)[1]) - Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  
  X_coeff <- matrix(runif(X_cols*n_cols, -1, 1), ncol = n_cols)
  Z_coeff <- matrix(runif(Z_cols*n_rows, -1, 1), ncol = n_rows)
  D_coeff <- matrix(runif(n_rows*10, -1, 1), nrow = n_rows) %*% matrix(runif(n_cols*10, -1, 1), ncol = n_cols)
  
  Y_data <- X_data%*%X_coeff + t(Z_data%*%Z_coeff) + Pxp %*% D_coeff %*% Pzp + E_data
  
  #*************
  output <- Y_data
  row_input <- X_data
  col_input <- Z_data
  mask <- W_data
  
  
  source("soft_mc.R")
  source("mao_original_x.R")
  source("mao_original_z.R")
  source("mao_coop.R")
  
  sink("sim1_07maoxz.txt", append = TRUE)
  print(c(round(mao_error_coop,4), round(mao_error_x,4), round(mao_error_z,4), round(soft_error,4)))
  
  sink()
  
}


