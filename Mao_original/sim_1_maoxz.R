
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/Mao")
rm(list = ls())
library(softImpute)
library(purrr)
source("SMC_functions6.R")
seed = 2023
set.seed(seed)
n_rows <- 300
n_cols <- 300
coll <- TRUE

X_cols <- 5
Z_cols <- 5

missing_prob <- 0.8
#W_data <- matrix( rbinom(n_rows*n_cols, 1, (1 - missing_prob) ) , nrow = n_rows)

### sim1 data
#**************

sink("sim1_07maoxz.txt", append = TRUE)
print(c("mao_error_coop", "mao_error_x", "mao_error_z", "soft_error"))
sink()
print("-----------------------------------------------------------------")
for(iter_i in 1:1){
  
  print(iter_i)
  # EDIT: Adding my simulation code to produce consistent results
  
  set.seed(seed=seed)
  X_data <- matrix(rnorm(n_rows*X_cols), ncol = X_cols)
  Z_data <- matrix(rnorm(n_cols*Z_cols), ncol = Z_cols)
  E_data <- matrix(rnorm(n_rows*n_cols), ncol = n_cols)
  # normalize X
  X_coeff <- matrix(runif(X_cols*n_cols, -1, 1), ncol = n_cols)
  Z_coeff <- matrix(runif(Z_cols*n_rows, -1, 1), ncol = n_rows)
  D_coeff <- matrix(runif(n_rows*Z_cols, -1, 1), nrow = n_rows) %*% matrix(runif(n_cols*Z_cols, -1, 1), ncol = n_cols)
  
  if(coll == TRUE){
    X_data[,2]  <- X_data[,1] + rnorm(n_rows, mean = 0, sd = 0.001)
    Z_data[,2]  <- Z_data[,1] + rnorm(n_cols, mean = 0, sd = 0.001)
  }
  Px <- X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  Pz <- Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  Pxp <- diag(dim(X_data)[1]) - X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  Pzp <- diag(dim(Z_data)[1]) - Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  
  W_data <- matrix( rbinom(n_rows*n_cols, 1, (1 - missing_prob) ) , nrow = n_rows)
  Y_data <- X_data%*%X_coeff + t(Z_data%*%Z_coeff) + Pxp %*% D_coeff %*% Pzp + E_data
  
  
  #A <- matrix(rnorm(n_rows*10), nrow = n_rows)
  #B <- matrix(rnorm(n_cols*10), ncol = n_cols)
  #C <- matrix(rnorm(n_rows*10), nrow = n_rows)
  #D <- matrix(rnorm(n_cols*10), ncol = n_cols)
  
  # END of EDIT
  
  
  # X_data <- matrix(rnorm(n_rows*X_cols), ncol = X_cols)
  # Z_data <- matrix(rnorm(n_cols*Z_cols), ncol = Z_cols)
  # E_data <- matrix(rnorm(n_rows*n_cols), ncol = n_cols)
  # A <- matrix(rnorm(n_rows*10), nrow = n_rows)
  # B <- matrix(rnorm(n_cols*10), ncol = n_cols)
  # C <- matrix(rnorm(n_rows*10), nrow = n_rows)
  # D <- matrix(rnorm(n_cols*10), ncol = n_cols)
  # Px <- X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  # Pz <- Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  # Pxp <- diag(dim(X_data)[1]) - X_data %*% solve(t(X_data) %*% X_data) %*% t(X_data)
  # Pzp <- diag(dim(Z_data)[1]) - Z_data %*% solve(t(Z_data) %*% Z_data) %*% t(Z_data)
  # 
  # X_coeff <- matrix(runif(X_cols*n_cols, -1, 1), ncol = n_cols)
  # Z_coeff <- matrix(runif(Z_cols*n_rows, -1, 1), ncol = n_rows)
  # D_coeff <- matrix(runif(n_rows*10, -1, 1), nrow = n_rows) %*% matrix(runif(n_cols*10, -1, 1), ncol = n_cols)
  # 
  # Y_data <- X_data%*%X_coeff + t(Z_data%*%Z_coeff) + Pxp %*% D_coeff %*% Pzp + E_data
  
  #*************
  output <- Y_data
  row_input <- X_data
  col_input <- Z_data
  mask <- W_data
  
  
  #source("soft_mc.R")
  #source("mao_original_x.R") 
  #source("mao_original_z.R")
  max_iter <- 10
  alpha <- 1  # the agreement penalty parameter
  source("mao_coop.R")
  alpha <- 0.9  # the agreement penalty parameter
  source("mao_coop.R")
  alpha <- 0.7  # the agreement penalty parameter
  source("mao_coop.R")
  alpha <- 0.5  # the agreement penalty parameter
  source("mao_coop.R")
  alpha <- 0.3  # the agreement penalty parameter
  source("mao_coop.R")
  print("end ----------------")
  
  sink("sim1_07maoxz.txt", append = TRUE)
  #print(c(round(mao_error_coop,4), round(mao_error_x,4), round(mao_error_z,4), round(soft_error,4)))
  sink()
  
}


