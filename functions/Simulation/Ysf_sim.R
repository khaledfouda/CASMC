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
            informative_cov_prop = .5) {
      if (!is.null(seed))
         set.seed(seed = seed)
      X <- matrix(rnorm(n1 * m1), ncol = m1)
      Z <- matrix(rnorm(n2 * m2), ncol = m2)
      E <- matrix(rnorm(n1 * n2), ncol = n2)
      # normalize X
      #X <- normalize_matrix(X)
      beta.x <- matrix(runif(m1 * n2), ncol = n2)
      beta.z <- matrix(runif(m2 * n1), ncol = n1)
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
         if (cov_eff) {
            O <- (X %*% beta.x) + P_bar_X %*% M
         } else{
            ncov_to_keep = round(informative_cov_prop * m1)
            if (ncov_to_keep < m1) {
               beta.x[(ncov_to_keep + 1):m1, ] <- 0
               # beta.z[(ncov_to_keep+1):m2,] <- 0
               xbeta = X[, 1:ncov_to_keep] %*% beta.x[1:ncov_to_keep, ]
            } else{
               beta.x <- matrix(0, m1, ncol = n2)
               xbeta <- X %*% beta.x
            }
            beta.z <- matrix(0, m2, ncol = n1)
            O <- xbeta + P_bar_X %*% M
         }
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
