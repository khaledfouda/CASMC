CASMC3_kfold_M_3 <-
 function(Y,
          # y_train is expected to be Incomplete
          X,
          obs_mask,
          n_folds = 3,
          # y: a final full-fit if provided. Expected to be Incomplete
          lambda.beta = .Machine$double.eps,
          learning.rate = 0.001,
          beta.iter.max = 20,
          # provide this if you need L2 regularization.
          error_function = error_metric$rmse,
          # tuning parameters for lambda
          lambda.factor = 1 / 4,
          lambda.init = NULL,
          n.lambda = 40,
          # tuning parameters for J
          rank.init = 2,
          rank.limit = 30,
          rank.step = 1,
          pct = 0.98,
          # laplacian parameters
          lambda.a = 0,
          S.a = NULL,
          lambda.b = 0,
          S.b = NULL,
          # stopping criteria
          early.stopping = 3,
          thresh = 1e-6,
          maxit = 100,
          # trace parameters
          trace = FALSE,
          print.best = TRUE,
          quiet = FALSE,
          # initial values.
          warm = NULL,
          # seed
          seed = NULL) {
  if (!is.null(seed))
   set.seed(seed)
  
  # prepare the sequence of lambda (nuclear regularization hyperparameter)
  if (is.null(lambda.init))
   lambda.init <-
    lambda0.cov_splr(Y, reduced_hat_decomp.H(X)) * lambda.factor
  lamseq <- seq(from = lambda.init,
                to = .Machine$double.eps,
                length = n.lambda)
  #-----------------------------------------------------------------------
  # prepare the folds
  folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, obs_mask)
  fold_data <- lapply(1:n_folds, function(i) {
   valid_mask = folds[[i]]
   valid_ind = valid_mask == 0 & obs_mask == 1
   
   Y_train = Y * valid_mask
   Y_train[valid_mask == 0] = NA
   Y_train = as(Y_train, "Incomplete")
   
   Y_valid = Y[valid_ind]
   
   valid_mask[valid_ind] = 1
   valid_mask[!valid_ind] = NA
   valid_mask <- as(valid_mask, "Incomplete")
   virow = valid_mask@i
   vpcol = valid_mask@p
   valid_mask <- NULL
   
   list(
    Y_train = Y_train,
    Y_valid = Y_valid,
    virow = virow,
    vpcol = vpcol
   )
  })
  #------------------------------------------------
  #-----------------------------------------------------------------------
  #---------------------------------------------------------------------
  warm <- NULL
   #-----------------------------------------------------------------------
   rank.max <- rank.init
   best_fit <- list(error = Inf)
   counter <- 0
   fold_counter <- 0
   #---------------------------------------------------------------------
   for (i in seq(along = lamseq)) {
    fold <- (fold_counter %% n_folds) + 1
    data = fold_data[[fold]]
    fold_counter = fold_counter + 1
    
    
    fiti <-
     CASMC3_fit(
      y = data$Y_train,
      X = X,
      J = rank.max,
      lambda.M = lamseq[i],
      learning.rate = learning.rate,
      beta.iter.max = beta.iter.max,
      lambda.beta = lambda.beta,
      lambda.a = lambda.a,
      S.a = S.a,
      lambda.b = lambda.b,
      S.b = S.b,
      warm.start = warm,
      trace.it = F,
      thresh = thresh,
      maxit = maxit
     )
    
    #--------------------------------------------------------------
    # predicting validation set and xbetas for next fit:
    XbetaValid = suvC(X, t(fiti$beta), data$virow, data$vpcol)
    MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), data$virow, data$vpcol)
    #--------------------------------------------
    err = error_function(MValid + XbetaValid, data$Y_valid)
    #rank <- sum(round(fiti$d, 4) > 0)
    # newly added, to be removed later
    var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
    cum_var = cumsum(var_explained)
    rank  <- which(cum_var >= pct)[1]
    
    warm <- fiti # warm start for next
    #---------------------------------------------------------------------
    if (trace == TRUE)
     print(sprintf(
      paste0(
       "%2d lambda.M = %.3f, rank.max = %d  ==>",
       " rank = %d, error = %.5f, niter/fit = %d [Beta(lambda=%.3f)]"
      ),
      i,
      lamseq[i],
      rank.max,
      rank,
      err,
      fiti$n_iter,
      lambda.beta
     ))
    #-------------------------
    # register best fir
    if (err < best_fit$error) {
     best_fit$error = err
     best_fit$rank_M = rank
     best_fit$lambda.M = lamseq[i]
     best_fit$rank.max = rank.max
     best_fit$fit = fiti
     best_fit$iter = i
     counter = 0
    } else
     counter = counter + 1
    if (counter >= early.stopping) {
     if (trace)
      print(
       sprintf(
        "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
        counter
       )
      )
     break
    }
    # compute rank.max for next iteration
    rank.max <- min(rank + rank.step, rank.limit)
   }
  #---------------------------------------------------------------------
  # fit one last time full model, if the train/valid is provided
  # Y[Y == 0] = NA
  # Y <- as(Y, "Incomplete")
  # best_fit$fit <-
  #  CASMC3_fit(
  #   y = Y,
  #   X = X,
  #   J = best_fit$rank.max,
  #   lambda.M = best_fit$lambda.M,
  #   learning.rate = learning.rate,
  #   beta.iter.max = beta.iter.max,
  #   lambda.beta = lambda.beta,
  #   lambda.a = lambda.a,
  #   S.a = S.a,
  #   lambda.b = lambda.b,
  #   S.b = S.b,
  #   warm.start = best_fit$fit,
  #   trace.it = F,
  #   thresh = thresh,
  #   maxit = maxit
  #  )
  # 
  
  best_fit$lambda.beta = lambda.beta
  best_fit$lambda.a = lambda.a
  best_fit$lambda.b = lambda.b
  best_fit$learning.rate = learning.rate
  best_fit$rank_M = best_fit$fit$J
  return(best_fit)
 }

#---------------------------------------------------------------------------------
