#----------------------------------------------------------------------
# This functions generates either model M1 or M2 from Yousef's simulation. Set model=1 for M1 or model=2 for M2.
utils$generate_simulated_data <- 
   generate_simulated_data <-
   function(n = 300,
            m = 400,
            r = 10,
            k = 10,
            missing_prob = 0.8,
            coll = FALSE,
            half_discrete = FALSE,
            informative_cov_prop = 1,
            prepare_for_fitting = FALSE,
            mar_sparse = FALSE,
            mv_beta = TRUE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed = seed)
      #-----------------------------------
      X <- matrix(rnorm(n * k), ncol = k)
      if (half_discrete) {
         ndisc = ceiling(k / 2) 
         sampled_disc = sample(1:k, ndisc, FALSE)
         X[, sampled_disc] <- rbinom(n * ndisc, 1, 0.3)
      }
      # if collinearity is needed, make the correlation between the first two columns in X and Z between (.99,1)
      # the actual correlation value is very close to 0.999999%
      if (coll == TRUE) {
         cor(X[, 1], X[, 1] + rnorm(n, mean = 0, sd = 1.5))
         X[, 2]  <- X[, 1] + rnorm(n, mean = 0, sd = 1.5)
      }
      #---------------------------
      if (mv_beta) {
         beta_means <- runif(k, 1, 3) * sample(c(-1, 1), k, TRUE)
         beta_vars <- runif(k, 0.5, 1) ^ 2
         beta <-
            t(mvrnorm(m, beta_means, diag(beta_vars, k, k)))
      } else
         beta <-
            matrix(runif(k * m, 1, 2) * sample(c(1, 1), k * m, TRUE), nrow = k)
      U <- matrix(runif(n * r), ncol = r)
      V <- matrix(runif(r * m), nrow = r)
      P_X = X %*% solve(t(X) %*% X) %*% t(X)
      P_bar_X = diag(1, n, n) - P_X
      
      M <- P_bar_X %*% U %*% V
      
      W <- matrix(rbinom(n * m, 1, (1 - missing_prob)) , nrow = n)
      
      if (informative_cov_prop < 1) {
         if (mar_sparse && informative_cov_prop > 0) {
            indices_to_zero <- sample(1:length(beta), (1-informative_cov_prop) * length(beta))
            beta[indices_to_zero] <- 0
         } else{
            ncov_to_keep = round(informative_cov_prop * k)
            
            if (ncov_to_keep <= 0) {
               beta <- matrix(0, k, ncol = m)
            } else if (ncov_to_keep < k) {
               sampled_covars_to_remove <-
                  sample(1:k, k - ncov_to_keep, replace = FALSE)
               beta[sampled_covars_to_remove,] <- 0
            }
         }
      }
      
      O <- X %*% beta +  M
      # random noise to Y
      sqrt_signal <-
         sqrt(sum((O - mean(O)) ^ 2) / (n * m - 1)) # set SNR=1
      E <-
         matrix(rnorm(n * m, mean = 0, sd = 1), ncol = m)
      Y <- (O + E) * W
      rank <- qr(O)$rank
      
      fit_data <- NULL
      if (prepare_for_fitting) {
         fit_data <- list()
         fit_data$W_valid <- matrix.split.train.test(W, testp = 0.2)
         train = (Y * fit_data$W_valid)
         fit_data$train = to_incomplete(train)
         fit_data$valid = Y[fit_data$W_valid == 0]
         fit_data$Y = to_incomplete(Y)
      }
      
      return(list(
         O = O,
         W = W,
         X = X,
         Y = Y,
         beta = beta,
         M = M,
         fit_data = fit_data,
         rank = rank
      ))
      
      
      #-----------------------------------------------------------------------------
      #---------------------------------------------------------------------
      
   }
