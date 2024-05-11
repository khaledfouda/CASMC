
#------------------------------------------------------------------------------------------
CASMC_cv_holdout_with_r <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            lambda.beta,
            lambda.M,
            J,
            y = NULL,
            r_min = 0,
            r_max = X_r$rank,
            error_function = error_metric$rmse,
            trace = FALSE,
            print.best = TRUE,
            thresh = 1e-6,
            maxit = 100,
            lambda.a = 0,
            S.a = NULL,
            lambda.b = 0,
            S.b = NULL,
            warm = NULL,
            track_r = FALSE,
            pct = 0.98,
            quiet = FALSE,
            seed = NULL) {
      r_seq <- (max(r_min, 0)):(min(X_r$rank, r_max))
      Xterms = GetXterms(X_r$X, lambda.beta)
      best_score = Inf
      best_fit = NULL
      W_valid[W_valid == 1] = NA
      W_valid[W_valid == 0] =  1
      W_valid <- as(W_valid, "Incomplete")
      virow = W_valid@i
      vpcol = W_valid@p
      W_valid = NULL
      
      results <- list()
      for(r in r_seq){
      
         fiti <-
            CASMC_fit_rank(
               y = y_train,
               X = X_r$X,
               Xterms = Xterms,
               r = r,
               J = J,
               lambda = lambda.M,
               warm.start = warm,
               trace.it = F,
               thresh = thresh,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               final.svd = T,
               maxit = maxit
            )
         
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         #XbetaValid = suvC(X_r$X , t(fiti$beta), virow, vpcol)
         Beta = fiti$Beta
         XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), virow, vpcol)
         MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
         err = error_function(MValid + XbetaValid, y_valid)
         #--------------------------------------------
         if(err < best_score){
            best_score = err
            best_fit = fiti
            best_r = r
         }
         if(track_r) print(paste("r = ",r, " - error = ", err))
         
      }
      
      # fit one last time full model, if the train/valid is provided
      if (!is.null(y)) {
         stopifnot(inherits(y, "dgCMatrix"))
         best_fit <-
            CASMC_fit_rank(
               y = y,
               X = X_r$X,
               Xterms = Xterms,
               r = r,
               J = J,
               lambda = lambda.M,
               warm.start = best_fit,
               thresh = thresh,
               maxit = maxit,
               trace.it = F,
               lambda.a = lambda.a,
               S.a = S.a,
               lambda.b = lambda.b,
               S.b = S.b,
               final.svd = T
            )
         
      }
      
      
      list(
         best_score = best_score,
         best_fit = best_fit,
         best_r = best_r,
         lambda.beta = lambda.beta,
         lambda.M = lambda.M,
         J=J
      )
      
      
   }


