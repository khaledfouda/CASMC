
# Make sure to run Mao_import_lib.R before running this file.

Mao_Coop_fit <- function(A, X, Z, W, maxiter=100, epsilon=1e-6,
                         rho = 1, # the agreement penalty
                         n_folds=5,
                         lambda.1_grid = seq(0,2,length=10),
                         lambda.2_grid = seq(.9, 0, length=10),
                         alpha_grid = seq(0.992, 1, length=5),
                         numCores = 4,
                         n1n2_optimized = TRUE,
                         theta_estimator = theta_default,
                         seed = 2023){
   #' Input:
   #'       A: The response matrix of dimension n1 by n2. We assume that Y = A * W
   #'       X: The row covariates of dimension n1 by m1
   #'       Z: The column covariates of dimension n2 by m2
   #'       W: Is the mask matrix with the same dimension as Y
   #'       niter:
   #' 
   #' Output:
   #'       A.hat
   #' ----------------------------------------------------------------
   Y = A * W
   Y_tr = t(Y)
   A_tr = t(A)
   W_tr = t(W)
   # to obtain initial estimates, we run the MAO function on each covariate matrix individually
   best_param = Mao.cv(A, X, Y, W, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                         numCores, n1n2_optimized, theta_estimator)$best_parameters
   A.hat_x = Mao.fit(Y, X, W, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                     n1n2_optimized, theta_estimator )$A_hat
   
   
   best_param = Mao.cv(A_tr, Z, Y_tr, W_tr, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                         numCores, n1n2_optimized, theta_estimator)$best_parameters
   A.hat_z = t(Mao.fit(Y_tr, Z, W_tr, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                       n1n2_optimized, theta_estimator )$A_hat)
   
   # reporting test errors
   test_error_x = test_error(A.hat_x[W==0], A[W==0])
   test_error_z = test_error(A.hat_z[W==0], A[W==0])
   #print(test_error_z)
   test_error_avg = test_error(((A.hat_x+A.hat_z)/2)[W==0], A[W==0])
   print(paste("Test error using only X is:",
         round(test_error_x,5),
         ", and using only Z:",
         round(test_error_z,5),
         ", and using their average:",
         round(test_error_avg,5)
   ))
   
   
   # We assume to row covariates are minimized first for now.
   iter = 0
   n1 = dim(Y)[1]
   n2 = dim(Y)[2]
   A.hat = A.hat_x + A.hat_z
   rho.1 = 1/(1+rho)
   rho.2 = (1-rho)
   diff = Inf
   
   while(iter < maxiter & diff > epsilon){
      
      # 1. row covariates, (Y1)   
      y.star = rho.1 * (Y - rho.2 * A.hat_z)
      best_param = Mao.cv(A, X, y.star, W, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                            numCores, n1n2_optimized, theta_estimator)$best_parameters
      A.hat_x = Mao.fit(y.star, X, W, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                        n1n2_optimized, theta_estimator )$A_hat
      # 2. Column covariates (Y2)
      y.star = rho.1 * (Y_tr - rho.2 * t(A.hat_x))
      best_param = Mao.cv(A_tr, Z, y.star, W_tr, n_folds, lambda.1_grid, lambda.2_grid, alpha_grid, seed,
                            numCores, n1n2_optimized, theta_estimator)$best_parameters
      A.hat_z = t(Mao.fit(y.star, Z, W_tr, best_param$lambda.1, best_param$lambda.2, best_param$alpha, 
                          n1n2_optimized, theta_estimator )$A_hat)
      # update stopping criteria
      diff = sqrt(mean((A.hat - A.hat_x - A.hat_z)**2))
      A.hat = A.hat_x + A.hat_z
      iter = iter + 1
      test_error_Coop = test_error(A.hat[W==0], A[W==0])
      print(paste("Iteration",iter, "- test MSE =",test_error_Coop, "- diff =",diff))
   }
   if(iter >= maxiter){
      print("Reached Max iterations before converging.")
   }else
      print("Converged.")
   
   test_error_Coop = test_error(A.hat[W==0], A[W==0])
   print(paste("Test error resulted from applying the Coop Model is",test_error_Coop))
   return(A.hat)
} 
