#' Covariate-Adjusted-Sparse-Matrix-completion
#' Cross-Validation function with holdout method (training / validation)
#'
#' @param y_train A sparse matrix of class incomplete (training)
#' @param X_r A list containing SvdH and X (see fit function for more details)
#' @param y_valid A vector of the validation set.
#' @param W_valid An indicator matrix where 0 = validation and 1 = train/test
#' @param y (optional) A sparse matrix of class Incomplete of the train and validation combined. A final fit step will be
#'                            applied if y is provided.
#' @return A list of u,d,v of M, Beta, and a vector of the observed Xbeta
#' @examples
#'  CASMC_fit(y,X,J=5)
#' @export
#'
CASMC_cv_holdout <-
   function(y_train,
            X_r,
            y_valid,
            W_valid,
            y = NULL,
            error_function = RMSE_error,
            lambda.factor = 1 / 4,
            lambda.init = NULL,
            n.lambda = 20,
            trace = FALSE,
            print.best = TRUE,
            early.stopping = 1,
            thresh = 1e-6,
            maxit = 100,
            rank.init = 2,
            rank.limit = 30,
            rank.step = 2,
            warm = NULL,
            quiet = FALSE) {
      # prepare the sequence of lambda (nuclear regularization hyperparameter)
      if (is.null(lambda.init))
         lambda.init <-
            lambda0.cov_splr(y_train, X_r$svdH) * lambda.factor
      lamseq <- seq(from = lambda.init,
                    to = 0,
                    length = n.lambda)
      #----------------------------------------------------
      stopifnot(inherits(y_train, "dgCMatrix"))
      # we only need the indices for validation from W_valid
      W_valid[W_valid == 1] = NA
      W_valid[W_valid == 0] =  1
      W_valid <- as(W_valid, "Incomplete")
      virow = W_valid@i
      vpcol = W_valid@p
      W_valid = NULL
      #------------------------------------------------
      Xterms = GetXterms(X_r$X)
      #-----------------------------------------------------------------------
      rank.max <- rank.init
      best_fit <- list(error = Inf)
      counter <- 0
      #---------------------------------------------------------------------
      for (i in seq(along = lamseq)) {
         fiti <-
            CASMC_fit(
               y = y_train,
               X = X_r$X,
               svdH = X_r$svdH,
               Xterms = Xterms,
               J = rank.max,
               lambda = lamseq[i],
               warm.start = warm,
               trace = F,
               thresh = thresh,
               init = "naive",
               final.svd = T,
               maxit = maxit
            )
         
         #--------------------------------------------------------------
         # predicting validation set and xbetas for next fit:
         Beta = fast.svd(as.matrix(fiti$beta))
         XbetaValid = suvC(X_r$X %*% Beta$v, t(Beta$d * t(Beta$u)), virow, vpcol)
         MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
         #--------------------------------------------
         err = error_function(MValid + XbetaValid, y_valid)
         rank <- sum(round(fiti$d, 4) > 0)
         warm <- fiti # warm start for next
         #---------------------------------------------------------------------
         if (trace == TRUE)
            print(sprintf(
               paste0(
                  "%2d lambda=%9.5g, rank.max = %d  ==>",
                  " rank = %d, error = %.5f, niter/fit = %d"
               ),
               i,
               lamseq[i],
               rank.max,
               rank,
               err,
               fiti$n_iter
            ))
         #-------------------------
         # register best fir
         if (err < best_fit$error) {
            best_fit$error = err
            best_fit$rank_M = rank
            best_fit$lambda = lamseq[i]
            best_fit$rank.max = rank.max
            best_fit$fit = fiti
            best_fit$iter = i
            counter = 0
         } else
            counter = counter + 1
         if (counter >= early.stopping) {
            if (trace)
               print(sprintf(
                  "Early stopping. Reached Peak point. Performance didn't improve for the last %d iterations.",
                  counter
               ))
            break
         }
         # compute rank.max for next iteration
         rank.max <- min(rank + rank.step, rank.limit)
      }
      # fit one last time full model, if the train/valid is provided
      if (!is.null(y)) {
         stopifnot(inherits(y, "dgCMatrix"))
         best_fit$fit <-
            CASMC_fit(
               y = y,
               X = X_r$X,
               svdH = X_r$svdH,
               Xterms = Xterms,
               J = best_fit$rank.max,
               lambda = best_fit$lambda,
               warm.start = best_fit$fit1,
               thresh = thresh,
               maxit = maxit,
               trace = F,
               final.svd = T
            )
         
      }
      return(best_fit)
   }
#---
