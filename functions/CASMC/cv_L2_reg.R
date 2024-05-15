CASMC_cv_L2 <-
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
            # L2 parameters
            lambda.beta.grid = "default",
            track_beta = FALSE,
            max_cores = 8,
            # seed
            seed = NULL) {
      if (identical(lambda.beta.grid, "default"))
         lambda.beta.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *  seq(0, 10, length.out = 20)
      
      num_cores = length(lambda.beta.grid)
      if (length(lambda.beta.grid) > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(lambda.beta.grid) / 2))
      print(paste("Running on", num_cores, "cores."))
      
      results <- mclapply(lambda.beta.grid, function(lambda.beta) {
         Xterms = GetXterms(X, lambda.beta)
         fiti = CASMC_cv_nuclear(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = NULL,
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
         fiti$lambda.beta = lambda.beta
         fiti
      }, mc.cores = num_cores)
      
      
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (track_beta) {
         sapply(results, function(x)
            print(
               paste0(
                  "lambda.beta = ",
                  x$lambda.beta,
                  " - Val Err = ",
                  round(x$error, 5)
               )
            ))
      }
      if (print.best)
         print(
            paste(
               "Best fit: lambda_beta = ",
               best_fit$lambda.beta,
               " - Validation Error: ",
               best_fit$error
            )
         )
      
      return(best_fit)
      
      
   }
