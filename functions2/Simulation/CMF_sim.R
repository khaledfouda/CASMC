# This functions generates either model M1 or M2 from Yousef's simulation. Set model=1 for M1 or model=2 for M2.
generate_simulation_data_cmf <-
   function(n = 300,
            m = 300,
            k = 5,
            r = 10,
            missing_prob = 0.7,
            coll = FALSE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed = seed)
      X <- matrix(rnorm(n * k), ncol = k)
      B <- matrix(rnorm(m * r), ncol = r)
      C <- matrix(rnorm(r * k), ncol = k)
      # center X?
      #X <- normalize_matrix(X)
      O <- X %*% ginv(t(C)) %*% t(B)
      rank <- qr(O)$rank # = m + r
      #-----------------------------------------------------------------------------
         W <- matrix(rbinom(n * m, 1, (1 - missing_prob)) , nrow = n)
      #----------------------------------------------------------
      #----------------------------------------------------------------------
      # Computing epsilon as iid zero mean Gaussian with variance chosen such that the signal-to-noise ratio (SNR) is 1
      signal_O <- sum((O - mean(O)) ^ 2) / (n * m - 1)
      sigma_epsilon <- sqrt(signal_O)  # Since SNR = 1
      epsilon <-
         matrix(rnorm(n * m, mean = 0, sd = sigma_epsilon), n, m)
      #--------------------------------------------------------------------------------------------
      # Y is a corrupted, partially-observed version of A. Yij=0 for missing data [maybe consider making it NA]
      Y <- (O + epsilon) * W
      #---------------------------------------------------------------------
      return(list(
         O = O,
         W = W,
         X = X,
         Y = Y,
         rank = rank
      ))
   }
