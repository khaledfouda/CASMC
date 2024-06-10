CASMC2_cv_beta <-
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
            lambda.beta.grid = "default1",
            track = FALSE,
            max_cores = 8,
            # seed
            seed = NULL) {
      
      
      if(identical(lambda.beta.grid, "default1"))
         lambda.beta.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *  seq(10, 0, length.out = 20)
      if(identical(lambda.beta.grid, "default2"))
         lambda.beta.grid = seq(propack.svd(naive_fit(y, X, TRUE),1)$d/4,
                                0,
                                length.out = 20)
      
      
      num_cores = length(lambda.beta.grid)
      if (length(lambda.beta.grid) > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(lambda.beta.grid) / 2))
      print(paste("Running on", num_cores, "cores."))

      results = list()
      rank.beta = 3
      rank.beta.max = reduced_hat_decomp(dat$X, pct = 0.999)$rank #qr(X)$rank
      counter = 0
      best_fit <- list(error = Inf)
      fiti <- NULL
      
      for(i in seq(along = lambda.beta.grid)) {
      # results <- mclapply(lambda.beta.grid, function(lambda.beta) {
         fiti = CASMC2_cv_M(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            r = rank.beta,
            lambda.beta = lambda.beta.grid[i],
            
            # error_function = error_function,
            # lambda.factor = lambda.factor,
            # lambda.init = lambda.init,
            # n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            # early.stopping = early.stopping,
            # thresh = thresh,
            # maxit = maxit,
            # rank.init = rank.init,
            # rank.limit = rank.limit,
            # rank.step = rank.step,
            # pct = pct,
            warm = NULL,#warm,
            # lambda.a = lambda.a,
            # S.a = S.a,
            # lambda.b = lambda.b,
            # S.b = S.b,
            quiet = quiet,
            seed = seed
         )
      #    fiti$lambda.beta = lambda.beta
      #    fiti
      # }, mc.cores = num_cores, mc.cleanup = TRUE)
         err = fiti$error
         fiti = fiti$fit
         # var_explained = diag(fiti$db) ^ 2 / sum(diag(fiti$db) ^ 2)
         # cum_var = cumsum(var_explained)
         # rank.r  <- which(cum_var >= pct)[1]
         rank.r <- sum(round(diag(fiti$db), 4) > 0)
         
         warm <- fiti # warm start for next
         
         if (trace == TRUE)
            print(sprintf(
               paste0(
                  "%2d lambda.beta=%9.5g, rank.max.beta = %d  ==>",
                  " rank.beta = %d, error = %.5f, niter/fit = %d"
               ),
               i,
               lambda.beta.grid[i],
               rank.beta,
               rank.r,
               err,
               fiti$n_iter
            ))
         #-------------------------
         # register best fir
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_beta = rank.r
            best_fit$lambda.beta = lambda.beta.grid[i]
            best_fit$rank.beta.max = rank.beta
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
         rank.beta <- min(rank.r + 1, rank.beta.max)
         
      }
      # best_fit <-
      #    results[[which.min(sapply(results, function(x)
      #       x$error))]]
      # 
      # if (track) {
      #    sapply(results, function(x)
      #       print(
      #          paste0(
      #             "lambda.beta = ",
      #             x$lambda.beta,
      #             " - Val Err = ",
      #             round(x$error, 5)
      #          )
      #       ))
      # }
      # if (print.best)
      #    print(
      #       paste(
      #          "Best fit: lambda_beta = ",
      #          best_fit$lambda.beta,
      #          " - Validation Error: ",
      #          best_fit$error
      #       )
      #    )
      
      return(best_fit)
      
      
   }
