normalize_matrix <- function(X) {
   # NOT USED FOR NOW
   # Normalize rows
   X_row_normalized <- t(apply(X, 1, function(row) {
      (row - mean(row)) / sd(row)
   }))
   
   # Normalize columns
   X_col_normalized <- apply(X_row_normalized, 2, function(col) {
      (col - mean(col)) / sd(col)
   })
   
   return(X_col_normalized)
}


# This function generates Mao's simulation as defined in (M0) above. Dimensions, missing probability are given as parameters.
generate_simulation_data_mao <-
   function(n1 = 400,
            n2 = 400,
            m = 20,
            r = 10,
            method = "MAR",
            seed = NULL,
            miss.prob = 0.8,
            cov_eff = TRUE,
            discrete = FALSE,
            informative_cov_prop = 1) {
      #' Input:
      #'      n1, n2: are the dimensions of the A, Y, and B matrices
      #'      m: number of covariates
      #'      r: Second dimension of U and V where  B = P_bar_X U V^T
      #'      method: If "MAR", missing at random method is employed to compute the missing probability and W
      #'      if "UNI" then a uniform method is employed with probability equal to miss.prob
      #'      seed: random seed
      if (!is.null(seed))
         set.seed(seed = seed)
      X <- matrix(rnorm(n1 * m), ncol = m) #%>% normalize_matrix()
      if(discrete){
         ndisc = ceiling(m/2) +1
         X[, ndisc:m] <- rbinom(n1*(m-ndisc+1),1,0.4)
      }
      
      beta <- matrix(rnorm(m * n2), ncol = n2)
      U <- matrix(rnorm(n1 * r), ncol = r)
      V <- matrix(rnorm(n2 * r), ncol = r)
      P_X = X %*% solve(t(X) %*% X) %*% t(X)
      P_bar_X = diag(1, n1) - P_X
      M = P_bar_X %*% U %*% t(V)
      if (cov_eff) {
         O <- X %*% beta + M
      } else{
         ncov_to_keep = round(informative_cov_prop * m)
         if (ncov_to_keep < m) {
            beta[(ncov_to_keep + 1):m, ] <- 0
            xbeta = X[, 1:ncov_to_keep] %*% beta[1:ncov_to_keep, ]
         } else{
            beta <- matrix(0, m, n2)
            xbeta <- X %*% beta
         }
         O <- xbeta + M
      }
      
      rank <- qr(O)$rank # = m + r
      #-----------------------------------------------------------------------------
      stopifnot(method %in% c("MAR", "UNI"))
      if (method == "MAR") {
         # Simulation Theta using missing-at-random MAR method
         gamma <-
            rbind(matrix(rnorm(1, -1.5, 0.1), 1), matrix(rnorm(3, 0.3, 0.1), 3))
         # only the first 3 columns of X are used
         # theta is the inverse of the probability of missingness following a logistic model
         inclusion_prob = 1 / (1 + exp(-(cbind(1, X[, 1:3]) %*% gamma))) #  all values should be close to 0.2
         #theta = 1 / inclusion_prob
         # Wij=1 if A is observed following the probabilities of missingness defined above.
         W <- matrix(rbinom(n1 * n2, 1, inclusion_prob) , nrow = n1)
         # sum(W==0)/(n1*n2) should be equal to 0.8
      } else{
         W <- matrix(rbinom(n1 * n2, 1, (1 - miss.prob)) , nrow = n1)
      }
      #----------------------------------------------------------
      # Does fully observed Y = A (ie,  without noise?)? In that case ignore the code below.
      #----------------------------------------------------------------------
      # Computing epsilon as iid zero mean Gaussian with variance chosen such that the signal-to-noise ratio (SNR) is 1
      signal_O <- sum((O - mean(O)) ^ 2) / (n1 * n2 - 1)
      sigma_epsilon <- sqrt(signal_O)  # Since SNR = 1
      epsilon <-
         matrix(rnorm(n1 * n2, mean = 0, sd = sigma_epsilon), n1, n2)
      #--------------------------------------------------------------------------------------------
      # Y is a corrupted, partially-observed version of A. Yij=0 for missing data [maybe consider making it NA]
      Y <- (O + epsilon) * W
      #---------------------------------------------------------------------
      return(list(
         O = O,
         W = W,
         X = X,
         Y = Y,
         beta = beta,
         M = M,
         #theta = theta,
         #gamma = gamma,
         #inclusion_prob = inclusion_prob,
         rank = rank
      ))
   }
