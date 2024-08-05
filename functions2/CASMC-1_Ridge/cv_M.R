#---------------------------------------------------------------------------------------------
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
CASMC1_cv_M <-
  function(y_train,
           # y_train is expected to be Incomplete
           X,
           y_valid,
           # y_valid is a vector
           W_valid,
           y = NULL,
           # y: a final full-fit if provided. Expected to be Incomplete
           r = NULL,
           # provide this if you need rank restriction. if not null, L2 reg will be ignored
           Xterms = NULL,
           svdH = NULL,
           hpar = NULL,
           # provide this if you need L2 regularization.
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
           # seed
           seed = NULL) {
    if (!is.null(seed))
      set.seed(seed)
    if(is.null(svdH)) 
      svdH <-  utils$reduced_hat_decomp.H(X)
    # prepare the sequence of lambda (nuclear regularization hyperparameter)
    if (is.null(hpar$M$lambda.init))
      hpar$M$lambda.init <-
        utils$lambdaM.max(y_train, svdH) * hpar$M$lambda.factor
    lamseq <- seq(from = hpar$M$lambda.init,
                  to = 0,
                  length = hpar$M$n.lambda)
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
    if (is.null(Xterms))
      Xterms = utils$GetXterms(X)
    #-----------------------------------------------------------------------
    rank.max <- hpar$M$rank.init
    best_fit <- list(error = Inf, r = r)
    counter <- 0
    #---------------------------------------------------------------------
    for (i in seq(along = lamseq)) {
      fiti <-
        CASMC1_fit(
          y = y_train,
          X = X,
          Xterms = Xterms,
          J = rank.max,
          lambda = lamseq[i],
          lambda.a = hpar$laplacian$lambda.a,
          S.a = hpar$laplacian$S.a,
          lambda.b = hpar$laplacian$lambda.b,
          S.b = hpar$laplacian$S.b,
          warm.start = warm,
          trace.it = F,
          thresh = thresh,
          final.svd = T,
          maxit = maxit
        )
      
      #--------------------------------------------------------------
      # predicting validation set and xbetas for next fit:
      XbetaValid = suvC(X , t(fiti$beta), virow, vpcol)
      MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
      #--------------------------------------------
      err = error_function(MValid + XbetaValid, y_valid)
      #rank <- sum(round(fiti$d, 4) > 0)
      # newly added, to be removed later
      var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
      cum_var = cumsum(var_explained)
      rank  <- which(cum_var >= hpar$M$pct)[1]
      
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
      if (counter >= hpar$M$early.stopping) {
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
      rank.max <- min(rank + hpar$M$rank.step, hpar$M$rank.limit)
    }
    # fit one last time full model, if the train/valid is provided
    if (!is.null(y)) {
      stopifnot(inherits(y, "dgCMatrix"))
      best_fit$fit <-
        CASMC1_fit(
          y = y,
          X = X,
          Xterms = Xterms,
          J = best_fit$rank.max,
          lambda = best_fit$lambda,
          lambda.a = hpar$laplacian$lambda.a,
          S.a = hpar$laplacian$S.a,
          lambda.b = hpar$laplacian$lambda.b,
          S.b = hpar$laplacian$S.b,
          warm.start = best_fit$fit,
          trace.it = F,
          thresh = thresh,
          final.svd = T,
          maxit = maxit
        )
      
    }
    return(best_fit)
  }
#---