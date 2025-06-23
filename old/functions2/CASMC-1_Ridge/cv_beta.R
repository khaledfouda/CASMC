CASMC_Ridge_hparams <-
   list(
      M = list(
         # tuning parameters for lambda
         lambda.factor = 1 / 4,
         lambda.init = NULL,
         n.lambda = 20,
         # tuning parameters for J
         rank.init = 2,
         rank.limit = 30,
         rank.step = 2,
         pct = 0.98,
         early.stopping = 1
      ),
      beta = list(
         n.lambda = 40,
         lambda.grid = "default"
      ),
      laplacian = list(
         # laplacian parameters
         lambda.a = 0,
         S.a = NULL,
         lambda.b = 0,
         S.b = NULL
      )
   )




CASMC1_cv <-
   CASMC_Ridge_cv <-
   function(y_train,
            # y_train is expected to be Incomplete
            X,
            y_valid,
            # y_valid is a vector
            W_valid,
            y = NULL,
            hpar = CASMC_Ridge_hparams,
            # y: a final full-fit if provided. Expected to be Incomplete
            error_function = utils$error_metric$rmse,
            # tuning parameters for lambda
            #lambda.factor = 1 / 4,
            #lambda.init = NULL,
            #n.lambda = 20,
            # tuning parameters for J
            #rank.init = 2,
            #rank.limit = 30,
            #rank.step = 2,
            #pct = 0.98,
            # laplacian parameters
            #lambda.a = 0,
            #S.a = NULL,
            #lambda.b = 0,
            #S.b = NULL,
            # stopping criteria
            #early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            # trace parameters
            trace = FALSE,
            print.best = TRUE,
            quiet = FALSE,
            # initial values.
            warm = NULL,
            # L2 parameters
            
            track = FALSE,
            max_cores = 8,
            # seed
            seed = NULL) {
      if (identical(hpar$beta$lambda.grid, "default"))
         hpar$beta$lambda.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *  
            seq(0, 10, length.out = hpar$beta$n.lambda)
      
      num_cores = length(hpar$beta$lambda.grid)
      if (length(hpar$beta$lambda.grid) > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(hpar$beta$lambda.grid) / 2))
      print(paste("Running on", num_cores, "cores."))
      
      results <- mclapply(hpar$beta$lambda.grid, function(lambda.beta) {
         Xterms = utils$GetXterms(X, lambda.beta)
         fiti = CASMC1_cv_M(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            Xterms = Xterms,
            r = NULL,
            hpar = hpar,
            error_function = error_function,
            #lambda.factor = lambda.factor,
            #lambda.init = lambda.init,
            #n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            #early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            #rank.init = rank.init,
            #rank.limit = rank.limit,
            #rank.step = rank.step,
            #pct = pct,
            warm = NULL,
            #lambda.a = lambda.a,
            #S.a = S.a,
            #lambda.b = lambda.b,
            #S.b = S.b,
            quiet = quiet,
            seed = seed
         )
         fiti$lambda.beta = lambda.beta
         fiti
      }, mc.cores = num_cores, mc.cleanup = TRUE)
      
      if (track) {
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
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      if (print.best)
         print(
            paste(
               "Best fit: lambda_beta = ",
               best_fit$lambda.beta,
               " - Validation Error: ",
               best_fit$error
            )
         )
      best_fit$init_param <- hpar
      return(best_fit)
      
      
   }
