CASMC3_cv_beta <-
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
            trace = 0,
            print.best = TRUE,
            quiet = FALSE,
            # initial values.
            warm = NULL,
            # L1 parameters
            lambda.beta.grid = "default1",
            learning.rate = .001,
            beta.iter.max = 20,
            track = FALSE,
            max_cores = 8,
            # seed
            seed = NULL) {
      if (identical(lambda.beta.grid, "default1")){
         nf <- naive_fit(y_train, X)
         resids <- y_train - nf$M - X %*% nf$beta
         resids[y_train==0] <- 0
         term <- max((nf$beta / learning.rate) - t(X) %*% resids) 
         lambda.beta.grid = seq(term, .Machine$double.eps, length.out=20)
         
         # lambda.beta.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *
         #    seq(10, .Machine$double.eps, length.out = 20)
      
      }else if (identical(lambda.beta.grid, "default2"))
         lambda.beta.grid = seq(propack.svd(naive_fit(y, X, TRUE), 1)$d / 4,
                                .Machine$double.eps,
                                length.out = 20)
      
      
      num_cores = length(lambda.beta.grid)
      if (length(lambda.beta.grid) > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(lambda.beta.grid) / 2))
      print(paste("Running on", num_cores, "cores."))
      
      best_fit <- list(error = Inf)
      results <- mclapply(lambda.beta.grid, function(lambda.beta) {
         fiti = tryCatch({
            CASMC3_cv_M(
               y_train = y_train,
               X = X,
               y_valid = y_valid,
               W_valid = W_valid,
               y = y,
               learning.rate = learning.rate,
               lambda.beta = lambda.beta,
               beta.iter.max = beta.iter.max,
               error_function = error_function,
               lambda.factor = lambda.factor,
               lambda.init = lambda.init,
               n.lambda = n.lambda,
               trace = ifelse(trace == 2, TRUE, FALSE),
               print.best = print.best,
               thresh = thresh,
               maxit = maxit,
               rank.init = rank.M.init,
               rank.limit = rank.M.limit,
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
         }, error = function(e)
            list(
               error_message = e,
               error = 999999,
               lambda.beta = lambda.beta
            ))
         fiti
      }, mc.cores = num_cores, mc.cleanup = TRUE)
      
      # showing errors, if any,
      sapply(results, function(x) {
         if (!is.null(x$error_message))
            print(
               paste(
                  "Error encountered at lambda.beta = ",
                  round(x$lambda.beta, 3),
                  "with the following error message: ",
                  x$error_message
               )
            )
      })
      
      # extract best fit
      best_fit <-
         results[[which.min(sapply(results, function(x)
            x$error))]]
      
      # print all fit output:
      if (trace > 0) {
         sapply(results, function(x)
            print(
               sprintf(
                  paste0(
                     "<< lambda.beta = %.3f, error = %.5f, niter/fit = %.0f, M = [%.0f,%.3f] >> "
                  ),
                  x$lambda.beta,
                  x$error,
                  x$fit$n_iter,
                  x$fit$J,
                  x$fit$lambda.M
               )
            ))
      }
      
      # print best if
      if (print.best)
         print(
            sprintf(
               paste0(
                  "<< Best fit >> lambda.beta = %.3f, error = %.5f, niter/fit = %.0f, M = [%.0f,%.3f] > "
               ),
               best_fit$lambda.beta,
               best_fit$error,
               best_fit$fit$n_iter,
               best_fit$fit$J,
               best_fit$fit$lambda.M
            )
         )
      #-------------------------
      
      best_fit$hparams = data.frame(
         lambda.beta = best_fit$lambda.beta,
         lambda.M = best_fit$fit$lambda.M,
         learning.rate = learning.rate,
         rank.M = best_fit$fit$J
      )
      
      return(best_fit)
      
      
   }
