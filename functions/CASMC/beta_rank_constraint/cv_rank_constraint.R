# CASMC_cv_holdout_with_r <-
#    function(y_train,
#             X_r,
#             y_valid,
#             W_valid,
#             r_min = 0,
#             y = NULL,
#             error_function = RMSE_error,
#             lambda.factor = 1 / 4,
#             lambda.init = NULL,
#             n.lambda = 20,
#             trace = FALSE,
#             print.best = TRUE,
#             early.stopping = 1,
#             thresh = 1e-6,
#             maxit = 100,
#             rank.init = 2,
#             rank.limit = 30,
#             rank.step = 2,
#             warm = NULL,
#             track_r = FALSE,
#             quiet = FALSE) {
#       r_seq <- (X_r$rank):(r_min)#(X_r$rank):(r_min)
#       Xterms = GetXterms(X_r$X)
#       best_score = Inf
#       best_fit = NULL
#       warm = NULL
#
#       for (r in r_seq) {
#          fiti <- CASMC_cv_holdout(
#             y_train = y_train,
#             X_r = X_r,
#             y_valid = y_valid,
#             W_valid = W_valid,
#             y = y,
#             Xterms = Xterms,
#             r = r,
#             error_function = error_function,
#             lambda.factor = lambda.factor,
#             lambda.init = lambda.init,
#             n.lambda = n.lambda,
#             trace = trace,
#             print.best = print.best,
#             early.stopping = early.stopping,
#             thresh = thresh,
#             maxit = maxit,
#             rank.init = rank.init,
#             rank.limit = rank.limit,
#             rank.step = rank.step,
#             warm = warm,
#             quiet = quiet
#          )
#          #warm = fiti$fit
#          if (fiti$error < best_score) {
#             best_score = fiti$error
#             best_fit = fiti
#          }
#          if (track_r)
#             print(paste(r, "-", fiti$error))
#       }
#       return(best_fit)
#
#    }

#------------------------------------------------------------------------------------------
CASMC_cv_rank <-
 function(y_train,
          # y_train is expected to be Incomplete
          X,
          y_valid,
          # y_valid is a vector
          W_valid,
          y = NULL,
          # y: a final full-fit if provided. Expected to be Incomplete
          r = NULL,
          # provide this if you need rank restriction. if not null, L2 reg will be ignored
          Xterms = NULL,
          # provide this if you need L2 regularization.
          error_function = error_metric$rmse,
          # tuning parameters for lambda
          lambda.factor = 1 / 4,
          lambda.init = NULL,
          n.lambda = 20,
          # tuning parameters for J
          rank.init = 2,
          rank.limit = 30,
          rank.step = 2,
          pct = 0.98,
          # laplacian parameters
          lambda.a = 0,
          S.a = NULL,
          lambda.b = 0,
          S.b = NULL,
          # stopping criteria
          early.stopping = 1,
          thresh = 1e-6,
          maxit = 100,
          # trace parameters
          trace = FALSE,
          print.best = TRUE,
          quiet = FALSE,
          # initial values.
          warm = NULL,
          # rank constraint parameters
          rank_x = qr(X)$rank,
          r_min = 0,
          r_max = rank_x,
          track_r = FALSE,
          max_cores = 8,
          # seed
          seed = NULL
          ) {
  r_seq <- (max(r_min, 0)):(min(rank_x, r_max))
  num_cores = min(max_cores, length(r_seq))
  print(paste("Running on", num_cores, "cores."))
  Xterms = GetXterms(X)
  best_score = Inf
  best_fit = NULL
  results <- mclapply(r_seq, function(r) {
   CASMC_cv_nuclear(
    y_train = y_train,
    X = X,
    y_valid = y_valid,
    W_valid = W_valid,
    y = y,
    Xterms = Xterms,
    r = r,
    error_function = error_function,
    lambda.factor = lambda.factor,
    lambda.init = lambda.init,
    n.lambda = n.lambda,
    trace = trace,
    print.best = print.best,
    early.stopping = early.stopping,
    thresh = thresh,
    maxit = maxit,
    rank.init = rank.init,
    rank.limit = rank.limit,
    rank.step = rank.step,
    pct = pct,
    warm = NULL,
    lambda.a = lambda.a,
    S.a = S.a,
    lambda.b = lambda.b,
    S.b = S.b,
    quiet = quiet,
    seed = seed
   )
  }, mc.cores = num_cores)
  
  
  best_fit <-
   results[[which.min(sapply(results, function(x)
    x$error))]]
  
  if (track_r) {
   sapply(results, function(x)
    print(paste0("r = ",x$r, " - Val Err = ", round(x$error,5))))
  }
  
  return(best_fit)
  
  
 }




