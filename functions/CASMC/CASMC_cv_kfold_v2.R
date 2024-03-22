CASMC_cv_kfold_v2 <-
   function(Y,
            X_r,
            W,
            n_folds = 5,
            lambda.factor = 1 / 4,
            lambda.init = NA,
            n.lambda = 30,
            trace = FALSE,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 3,
            rank.limit = 20,
            rank.step = 2,
            warm = NULL,
            tol = 1,
            print.best = TRUE,
            seed = NULL) {
      if (!is.null(seed))
         set.seed(seed)
      #----------------------------------------------------
      lam0 <-
         ifelse(is.na(lambda.init),
                lambda0.cov_splr(Y, X_r$svdH) * lambda.factor,
                lambda.init)
      lamseq <- seq(from = lam0,
                    to = 1,
                    length = n.lambda)
      #-------------
      Xterms = GetXterms(X_r$X)
      #---------------------------------------------------------
      #-----------------------------------------------------------------------
      # prepare the folds
      folds <- k_fold_cells(nrow(Y), ncol(Y), n_folds, W)
      fold_data <- lapply(1:n_folds, function(i) {
         W_fold = folds[[i]]
         valid_ind = W_fold == 0 & W == 1
         
         Y_train = Y * W_fold
         Y_train[W_fold == 0] = NA
         Y_train = as(Y_train, "Incomplete")
         
         Y_valid = Y[valid_ind]
         
         W_fold[valid_ind] = 1
         W_fold[!valid_ind] = NA
         W_fold <- as(W_fold, "Incomplete")
         virow = W_fold@i
         vpcol = W_fold@p
         W_fold <- NULL
         
         list(
            Y_train = Y_train,
            Y_valid = Y_valid,
            virow = virow,
            vpcol = vpcol
         )
      })
      #---------------------------------------------------------------------------
      #n <- dim(Y)
      #m <- n[2]
      #n <- n[1]
      Y[Y == 0] = NA
      Y <- as(Y, "Incomplete")
      #---------------------------------------------------------------------------
      rank.max <- rank.init
      fiti = NULL
      counter = 0
      best_fit = list(error = Inf)
      niter1 <- 0
      #---------------------------------------------------------------------
      for (i in seq(along = lamseq)) {
         fiti <-
            CASMC_fit_v3(
               y = Y,
               X = X_r$X,
               svdH = X_r$svdH,
               trace = F,
               J = rank.max,
               thresh = thresh,
               lambda = lamseq[i],
               init = "naive",
               final.svd = T,
               maxit = maxit,
               warm.start = fiti,
               Xterms = Xterms
            )
         # #-----------------------------------------------
         err <- rank <- 0
         #-----------------------------------------------------
         for (fold in 1:n_folds) {
            data = fold_data[[fold]]
            #-----
            fiti <-
               CASMC_fit_v3(
                  y = data$Y_train,
                  X = X_r$X,
                  svdH = X_r$svdH,
                  trace = F,
                  J = rank.max,
                  thresh = thresh,
                  lambda = lamseq[i],
                  init = "naive",
                  final.svd = T,
                  maxit = maxit,
                  warm.start = fiti,
                  Xterms = Xterms
               )
            #--------------------------------------------------------------
            # predicting validation set and xbetas for next fit:
            Beta = fast.svd(as.matrix(fiti$beta))
            XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), data$virow, data$vpcol)
            MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), data$virow, data$vpcol)
            #--------------------------------------------
            err = err + test_error(MValid + XbetaValid, data$Y_valid)
            #--------------------------------------------
            rank <-
               rank + sum(round(fiti$d, 4) > 0) # number of positive sing.values
            if (trace) {
               niter1 <- niter1 + fiti$n_iter
            }
            #-----------------------------------------------------------------------------------
         }
         err = err / n_folds
         rank = as.integer(rank / n_folds)
         #------------------------------------------------
         if (trace) {
            print(sprintf(
               paste0(
                  "%2d lambda=%9.5g, rank.max = %d  ==>",
                  " rank = %d, error = %.5f, niter = %d"
               ),
               i,
               lamseq[i],
               rank.max,
               rank,
               err,
               niter1
            ))
            niter1 <- 0
         }
         #-----------------------------------------------------------------------
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$lambda = lamseq[i]
            best_fit$rank.max = rank.max
            best_fit$rank = rank
            best_fit$fit = fiti
            best_fit$iter = i
            counter = 0
         } else
            counter = counter + 1
         if (counter >= tol) {
            if (trace | print.best)
               print(sprintf(
                  "Performance didn't improve for the last %d iterations.",
                  counter
               ))
            break
         }
         #-------------------------------------------------------------------
         rank.max <- min(rank + rank.step, rank.limit)
         #----------------------------------------------------------------
         
      }
      if (print.best)
         print(
            sprintf(
               "lambda=%9.5g, rank.max = %d, error = %.5f",
               best_fit$lambda,
               best_fit$rank.max,
               best_fit$error
            )
         )
      #
      #    # one last fit!
      best_fit$fit <-
         CASMC_fit_v3(
            y = Y,
            X = X_r$X,
            svdH = X_r$svdH,
            trace = F,
            J = best_fit$rank.max,
            thresh = thresh,
            lambda = best_fit$lambda,
            init = "naive",
            final.svd = T,
            maxit = maxit,
            warm.start = best_fit$fit,
            Xterms = Xterms
         )
      
      
      return(best_fit)
   }
