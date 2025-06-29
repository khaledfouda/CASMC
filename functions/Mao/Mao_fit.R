#' The following function, $Mao.fit()$, estimates $\hat\beta$ and $\hat B$ as above with fixed (given) hyperparameters.
#'  I will try to use the same notations as above to avoid confusion.


# EDIT: These functions return the 1 / (probability of inclusion) NOT missingness.
Mao_weights <- list(
   binomial = function(X, W, ...) {
      # using logistic regression as indicated in (a1)
      n1 = dim(W)[1]
      n2 = dim(W)[2]
      theta_hat = matrix(NA, n1, n2)
      for (j in 1:n2) {
         model_data = data.frame(cbind(W[, j], X))
         model_fit = glm(X1 ~ ., family = binomial(), data = model_data)
         theta_hat[, j] = 1 / predict(model_fit, type = "response")
      }
      theta_hat[is.na(theta_hat) | is.infinite(theta_hat)] <- 0
      return(theta_hat)
   },
   column_avg = function(W, ...) {
      # A theta estimation function that selects the proportion of missing data within the same column
      # using formula (a2)
      n1 = dim(W)[1]
      n2 = dim(W)[2]
      theta_hat = matrix(NA, n1, n2)
      for (j in 1:n2) {
         theta_hat[, j] = n1 / sum(W[, j] == 1)
      }
      theta_hat[is.na(theta_hat) | is.infinite(theta_hat)] <- 0
      return(theta_hat)
   },
   uniform = function(W, ...) {
      # A theta estimation function that selects the proportion of missing data in the matrix
      # using formula (a3)
      n1 = dim(W)[1]
      n2 = dim(W)[2]
      theta_hat = matrix((n1 * n2) / sum(W == 1) , n1, n2)
      return(theta_hat)
   }
   
)

Mao_fit <-
   function(Y,
            X,
            W,
            lambda_1,
            lambda_2,
            alpha,
            n1n2_optimized = TRUE,
            return_rank = TRUE,
            theta_estimator = Mao_weights$binomial) {
      #
      #' ----------------------------------------------
      #' Input: Y: corrupted, partially observed A, (Y is assumed to be the product of Y*W)
      #'         Important: set missing values in Y to 0
      #'        W: Wij=1 if Aij is observed (ie, Yij !=0), and 0 otherwise
      #'         X: covariate matrix
      #'         lambda_1, lambda_2, alpha: hyperparameters
      #' ----------------------------------------------
      #' output: list of  A, Beta_hat, B_hat
      #' ----------------------------------------------
      n1 = dim(Y)[1]
      n2 = dim(Y)[2]
      m  = dim(X)[2]
      #yobs = W==1
      # The following two lines are as shown in (c) and (d)
      X.X = t(X) %*% X
      P_X = X %*% ginv(X.X) %*% t(X)
      P_bar_X = diag(1, n1, n1) - P_X
      
      if (n1n2_optimized == TRUE) {
         # we define the factor that will be used later:
         n1n2 = svd(X.X)$d[1]
      } else{
         n1n2 = n1 * n2 / 2
      }
      
      # The following part estimates theta (missingness probabilities)
      theta_hat = theta_estimator(W = W, X = X)
      # the following is the product of W * theta_hat * Y
      W_theta_Y = Y * theta_hat # * W
      
      # beta hat as (8)
      beta_hat = ginv(X.X + diag(n1n2 * lambda_1, m, m)) %*% t(X) %*% W_theta_Y
      # SVD decomposition to be used in (b)
      svdd = svd(P_bar_X %*% W_theta_Y)
      if (n1n2_optimized == TRUE) {
         # evaluation of  (b)
         n1n2 = svdd$d[1]
      } else{
         n1n2 = n1 * n2 / 2
      }
      T_c_D = svdd$u %*% (pmax(svdd$d - alpha * n1n2 * lambda_2, 0) * t(svdd$v))
      # B hat as in (11)
      B_hat = T_c_D / (1 + 2 * (1 - alpha) * n1n2 * lambda_2) 
      # computing the rank of B [Copied from Mao's code; Don't understand how it works.]
      # EQUIVALENT to  qr(B_hat)$rank + m   or   qr(A_hat)$rank
      # B is a low rank matrix
      #rank = sum(pmax(svdd$d - n1n2 * lambda_2 * alpha, 0) > 0) + m
      
      # Estimate the matrix as given in the model at the top
      A_hat = X %*% beta_hat + B_hat
      #A_hat[yobs] <- Y[yobs]
      rank = NULL
      if (return_rank)
         rank = qr(A_hat)$rank
      
      return(list(
         estimates = A_hat,
         M = B_hat,
         beta = beta_hat,
         rank = rank
      ))
   }
