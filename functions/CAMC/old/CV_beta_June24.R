



 CAMC_Lasso_cv <-
   function(y_train,
            # y_train is expected to be Incomplete
            X,
            y_valid,
            # y_valid is a vector
            W_valid,
            y = NULL,
            # y: a final full-fit if provided. Expected to be Incomplete
            hpar = CAMC_Lasso_hparams,
            
            # stopping criteria
            error_function = utils$error_metric$rmse,
            thresh = 1e-6,
            maxit = 100,
            # trace parameters
            trace = 0,
            print.best = TRUE,
            quiet = FALSE,
            # initial values.
            warm = NULL,
            
            track = FALSE,
            max_cores = 8,
            # seed
            seed = NULL) {
      
      if(identical(hpar$beta$learning.rate, "default"))
         hpar$beta$learning.rate <-  1 / sqrt(sum((t(X) %*% X)^2))
      
      if(is.null(hpar$beta$lambda.max)){
         
         nf <- naive_fit(y_train, X)
         resids <- y_train - nf$M - X %*% nf$beta
         resids[y_train==0] <- 0
         hpar$beta$lambda.max <- max((nf$beta / hpar$beta$learning.rate) - t(X) %*% resids) 
      }
         lambda.beta.grid = seq(hpar$beta$lambda.max, .Machine$double.eps, length.out=hpar$beta$n.lambda)
         
         
         # lambda.beta.grid = sqrt((ncol(y_train) * ncol(X)) / (nrow(y_train))) *
         #    seq(10, .Machine$double.eps, length.out = 20)
      
      
      
      num_cores = length(lambda.beta.grid)
      if (length(lambda.beta.grid) > max_cores)
         num_cores <-
         min(max_cores, ceiling(length(lambda.beta.grid) / 2))
      print(paste("Running on", num_cores, "cores."))
      
      best_fit <- list(error = Inf)
      results <- mclapply(lambda.beta.grid, function(lambda.beta) {
         fiti = tryCatch({
            CAMC3_cv_M(
               y_train = y_train,
               X = X,
               y_valid = y_valid,
               W_valid = W_valid,
               y = y,
               hpar = hpar,
               lambda.beta = lambda.beta,
               error_function = error_function,
               trace = ifelse(trace == 2, TRUE, FALSE),
               print.best = print.best,
               thresh = thresh,
               maxit = maxit,
               warm = NULL,
               
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
         learning.rate = hpar$beta$learning.rate,
         rank.M = best_fit$fit$J
      )
      best_fit$init_params = hpar
      return(best_fit)
      
      
   }
