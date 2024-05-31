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
           track = FALSE,
           max_cores = 8,
           # seed
           seed = NULL) {
    if(!is.null(seed)) set.seed(seed)
    r_seq <- (max(r_min, 0)):(min(rank_x, r_max))
    num_cores = min(max_cores, length(r_seq))
    print(paste("Running on", num_cores, "cores."))
    Xterms = GetXterms(X)
    svdH = reduced_hat_decomp.H(X)
    results <- mclapply(r_seq, function(r) {
    fiti <- tryCatch({
      CASMC_cv_nuclear(
        y_train = y_train,
        X = X,
        y_valid = y_valid,
        W_valid = W_valid,
        y = y,
        Xterms = Xterms,
        svdH = svdH,
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
    }, error = function(e) list(em=e, error=999999, r=r))
    fiti
    }, mc.cores = num_cores, mc.cleanup = TRUE)

    sapply(results, function(x) if(!is.null(x$em)) print(x$em))
    best_fit <-
      results[[which.min(sapply(results, function(x)
        x$error))]]
    if (track) {
      sapply(results, function(x)
        print(paste0(
          "r = ", x$r, " - Val Err = ", round(x$error, 5)
        )))
    }
    
    if (print.best)
      print(paste(
        "Best fit: r = ",
        best_fit$r,
        " - Validation Error: ",
        best_fit$error
      ))
    
    return(best_fit)
    
    
  }
