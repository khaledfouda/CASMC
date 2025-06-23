CASMC2_cv <-
   function(y_train,
            # y_train is expected to be Incomplete
            X,
            y_valid,
            # y_valid is a vector
            W_valid,
            y = NULL,
            rank.beta.init = 1,
            rank.beta.step = 1,
            rank.beta.limit = qr(X)$rank,
            lambda.beta.length = 20,
            lambda.beta.max = NULL,
            # y: a final full-fit if provided. Expected to be Incomplete
            error_function = error_metric$rmse,
            # tuning parameters for lambda
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            # tuning parameters for J
            rank.M.init = 2,
            rank.M.limit = 30,
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
      
      if(is.null(lambda.beta.max)){
         
      
      if (identical(lambda.beta.grid, "default1"))
         lambda.beta.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *
            seq(10, .Machine$double.eps, length.out = lambda.beta.length)
      
      
      if (identical(lambda.beta.grid, "default2"))
         lambda.beta.grid = seq(propack.svd(naive_fit(y, X, TRUE), 1)$d / 4,
                                .Machine$double.eps,
                                length.out = lambda.beta.length)
      }else
         lambda.beta.grid = seq(lambda.beta.max, 0, length.out=lambda.beta.length)
      
      # num_cores = length(lambda.beta.grid)
      # if (length(lambda.beta.grid) > max_cores)
      #    num_cores <-
      #    min(max_cores, ceiling(length(lambda.beta.grid) / 2))
      # print(paste("Running on", num_cores, "cores."))
      
      rank.max = rank.beta.init
      counter = 0
      best_fit <- list(error = Inf)
      
      for (i in seq(along = lambda.beta.grid)) {
         fiti = CASMC2_cv_M(
            y_train = y_train,
            X = X,
            y_valid = y_valid,
            W_valid = W_valid,
            y = y,
            r = rank.max,
            lambda.beta = lambda.beta.grid[i],
            
            error_function = error_function,
            lambda.factor = lambda.factor,
            lambda.init = lambda.init,
            n.lambda = n.lambda,
            trace = trace,
            print.best = print.best,
            #early.stopping = early.stopping,
            thresh = thresh,
            maxit = maxit,
            rank.init = rank.M.init,
            rank.limit = rank.M.limit,
            rank.step = rank.step,
            pct = pct,
            warm = warm,
            #warm,
            lambda.a = lambda.a,
            S.a = S.a,
            lambda.b = lambda.b,
            S.b = S.b,
            quiet = quiet,
            seed = seed
         )
         err = fiti$error
         fiti = fiti$fit
         rank <- sum(round(fiti$db, 6) > 0)
         # var_explained = diag(fiti$db) ^ 2 / sum(diag(fiti$db) ^ 2)
         # cum_var = cumsum(var_explained)
         # rank  <- which(cum_var >= pct)[1]
         
         warm <- fiti # warm start for next
         
         if (track == TRUE)
            print(sprintf(
               paste0(
                  "%2d lambda.beta=%9.5g, rank.max.beta = %d  ==>",
                  " rank.max = %d, error = %.5f, niter/fit = %d"
               ),
               i,
               lambda.beta.grid[i],
               rank.max,
               rank,
               err,
               fiti$n_iter
            ))
         #-------------------------
         # register best fit
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_beta = rank
            best_fit$lambda.beta = lambda.beta.grid[i]
            best_fit$rank.beta.max = rank.max
            best_fit$fit = fiti
            best_fit$iter = i
            counter = 0
         } else
            counter = counter + 1
         if (counter >= early.stopping) {
            if (track)
               print(
                  sprintf(
                     "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
                     counter
                  )
               )
            break
         }
         # compute rank.max for next iteration
         rank.max <- min(rank + rank.beta.step, rank.beta.limit)
         
      }
      best_fit$hparams = data.frame(lambda.beta = best_fit$lambda.beta,
                              lambda.M = best_fit$fit$lambda.M,
                              rank.beta = best_fit$rank_beta,
                              rank.M = best_fit$fit$J)
      best_fit$fit$beta = list(u = best_fit$fit$ub,
                               v = best_fit$fit$vb,
                               d = best_fit$fit$db^2)
      
      return(best_fit)
      
      
   }
