# This functions generates either model M1 or M2 from Yousef's simulation. Set model=1 for M1 or model=2 for M2.
generate_simulation_data_ysf <-
   function(model = 1,
            n1 = 300,
            n2 = 300,
            m1 = 5,
            m2 = 10,
            missing_prob = 0.7,
            coll = FALSE,
            seed = NULL,
            cov_eff = TRUE,
            discrete = FALSE,
            informative_cov_prop = 1) {
      if (!is.null(seed))
         set.seed(seed = seed)
      X <- matrix(rnorm(n1 * m1), ncol = m1)
      if (discrete) {
         ndisc = ceiling(m1 / 2) + 1
         X[, ndisc:m1] <- rbinom(n1 * (m1 - ndisc + 1), 1, 0.4)
      }
      Z <- matrix(rnorm(n2 * m2), ncol = m2)
      E <- matrix(rnorm(n1 * n2), ncol = n2)
      # normalize X
      #X <- normalize_matrix(X)
      beta.x <- matrix(runif(m1 * n2, 0.1, 1), ncol = n2)
      beta.z <- matrix(runif(m2 * n1, 0.1, 1), ncol = n1)
      M.x <- matrix(runif(n1 * m2), ncol = m2)
      M.z <- matrix(runif(m2 * n2), ncol = n2)
      M <- M.x %*% M.z
      
      # if collinearity is needed, make the correlation between the first two columns in X and Z between (.99,1)
      # the actual correlation value is very close to 95%
      if (coll == TRUE) {
         X[, 2]  <- X[, 1] + rnorm(n1, mean = 0, sd = 0.001)
         Z[, 2]  <- Z[, 1] + rnorm(n2, mean = 0, sd = 0.001)
      }
      
      
      P_X = X %*% solve(t(X) %*% X) %*% t(X)
      P_bar_X = diag(1, n1) - P_X
      P_Z = Z %*% solve(t(Z) %*% Z) %*% t(Z)
      P_bar_Z = diag(1, n2) - P_Z
      
      #X <- matrix(rnorm(n1 * m1, 5, 3), ncol = m1)
      
      W <-
         matrix(rbinom(n1 * n2, 1, (1 - missing_prob)) , nrow = n1)
      
      if (model == 1) {
         if (cov_eff) {
            O <- (X %*% beta.x) + t(Z %*% beta.z) + P_bar_X %*% M %*% P_bar_Z
         } else{
            beta.x <- matrix(0, m1, ncol = n2)
            beta.z <- matrix(0, m2, ncol = n1)
            O <-  P_bar_X %*% M %*% P_bar_Z
         }
         Y <- (O + E) * W
         rank <- qr(O)$rank
         return(list(
            O = O,
            W = W,
            X = X,
            Z = Z,
            Y = Y,
            beta.x = beta.x,
            beta.z = beta.z,
            M = M,
            rank = rank
         ))
      } else if (model == 2) {
         M <- P_bar_X %*% M
         if (cov_eff) {
            O <- (X %*% beta.x) + M
         } else{
            ncov_to_keep = round(informative_cov_prop * m1)
            if (ncov_to_keep < m1) {
               beta.x[(ncov_to_keep + 1):m1,] <- 0
               xbeta = X[, 1:ncov_to_keep] %*% beta.x[1:ncov_to_keep,]
            } else{
               beta.x <- matrix(0, m1, ncol = n2)
               xbeta <- X %*% beta.x
            }
            beta.z <- matrix(0, m2, ncol = n1)
            O <- xbeta +  M
         }
         #O = (O - mean(O)) / sd(O)
         Y <- (O + E) * W
         rank <- qr(O)$rank
         return(list(
            O = O,
            W = W,
            X = X,
            Y = Y,
            beta = beta.x,
            M = M,
            rank = rank
         ))
      }
      else{
         stop("Error: Unrecognized model.")
      }
      #-----------------------------------------------------------------------------
      #---------------------------------------------------------------------
      
   }
#----------------------------------------------------------------------
# This functions generates either model M1 or M2 from Yousef's simulation. Set model=1 for M1 or model=2 for M2.
generate_simulation_rows <-
   function(n = 300,
            m = 400,
            r = 10,
            k = 10,
            missing_prob = 0.8,
            coll = FALSE,
            half_discrete = FALSE,
            informative_cov_prop = 1,
            prepare_for_fitting = FALSE,
            mv_beta = TRUE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed = seed)
      #-----------------------------------
      X <- matrix(rnorm(n * k), ncol = k)
      if (half_discrete) {
         ndisc = ceiling(k / 2) + 1
         sampled_disc = sample(1:k, ndisc, FALSE)
         X[, sampled_disc] <- rbinom(n * (k - ndisc + 1), 1, 0.3)
      }
      # if collinearity is needed, make the correlation between the first two columns in X and Z between (.99,1)
      # the actual correlation value is very close to 0.999999%
      if (coll == TRUE) {
         X[, 2]  <- X[, 1] + rnorm(n, mean = 0, sd = 0.8)
         cor(X[,1], X[,1] + rnorm(n, mean = 0, sd = 0.8))
      }
      #---------------------------
      if (mv_beta) {
         beta_means <- runif(k, 1, 3) * sample(c(-1, 1), k, TRUE)
         beta_vars <- runif(k, 0.5, 1) ^ 2
         beta <-
            t(mvrnorm(m, beta_means, diag(beta_vars, k, k)))
      } else
         beta <- matrix(runif(k * m, 1, 2)*sample(c(1,1),k*m,TRUE), nrow = k)
      U <- matrix(runif(n * r), ncol = r)
      V <- matrix(runif(r * m), nrow = r)
      P_X = X %*% solve(t(X) %*% X) %*% t(X)
      P_bar_X = diag(1, n, n) - P_X
      
      M <- P_bar_X %*% U %*% V
      
      W <- matrix(rbinom(n * m, 1, (1 - missing_prob)) , nrow = n)
      
      ncov_to_keep = round(informative_cov_prop * k)
      
      if (ncov_to_keep <= 0) {
         beta <- matrix(0, k, ncol = m)
      } else if (ncov_to_keep < k) {
         sampled_covars_to_remove <-
            sample(1:k, k - ncov_to_keep, replace = FALSE)
         beta[sampled_covars_to_remove,] <- 0
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
