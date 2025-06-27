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
#'  CAMC_fit(y,X,J=5)
#' @export
#'
CAMC3_cv_M <-
  function(y_train,
           # y_train is expected to be Incomplete
           X,
           y_valid,
           # y_valid is a vector
           W_valid,
           y = NULL,
           # y: a final full-fit if provided. Expected to be Incomplete
           lambda.beta = .Machine$double.eps,
           hpar = hpar,
           error_function = utils$error_metric$rmse,
           # stopping criteria
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
    
    # prepare the sequence of lambda (nuclear regularization hyperparameter)
    if (is.null(hpar$M$lambda.init))
      hpar$M$lambda.init <-
        utils$lambdaM.max(y_train, utils$reduced_hat_decomp.H(X)) * hpar$M$lambda.factor
    lamseq <- seq(from = hpar$M$lambda.init,
                  to = .Machine$double.eps,
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
    #-----------------------------------------------------------------------
    rank.max <- hpar$M$rank.init
    best_fit <- list(error = Inf)
    counter <- 0
    #---------------------------------------------------------------------
    for (i in seq(along = lamseq)) {
      fiti <-
        CAMC3_fit(
          y = y_train,
          X = X,
          J = rank.max,
          lambda.M = lamseq[i],
          learning.rate = hpar$beta$learning.rate,
          beta.iter.max = hpar$beta$prox.iter.max,
          lambda.beta = lambda.beta,
          lambda.a = hpar$laplacian$lambda.a,
          S.a = hpar$laplacian$S.a,
          lambda.b = hpar$laplacian$lambda.b,
          S.b = hpar$laplacian$S.b,
          warm.start = warm,
          trace.it = F,
          thresh = thresh,
          maxit = maxit
        )
      
      #--------------------------------------------------------------
      # predicting validation set and xbetas for next fit:
      XbetaValid = suvC(X, t(fiti$beta), virow, vpcol)
      MValid = suvC(fiti$u, t(fiti$d * t(fiti$v)), virow, vpcol)
      #--------------------------------------------
      err = error_function(MValid + XbetaValid, y_valid)
      rank <- sum(round(fiti$d, 4) > 0)
      # newly added, to be removed later
      #var_explained = fiti$d ^ 2 / sum(fiti$d ^ 2)
      #cum_var = cumsum(var_explained)
      #rank  <- which(cum_var >= hpar$M$pct)[1]
      
      warm <- fiti # warm start for next
      #---------------------------------------------------------------------
      if (trace == TRUE)
        print(sprintf(
          paste0(
            "%2d lambda.M = %.3f, rank.max = %d  ==>",
            " rank = %d, error = %.5f, niter/fit = %d [Beta(lambda=%.3f)]"
          ),
          i,
          lamseq[i],
          rank.max,
          rank,
          err,
          fiti$n_iter,
          lambda.beta
        ))
      #-------------------------
      # register best fir
      if (err < best_fit$error) {
        best_fit$error = err
        best_fit$rank_M = rank
        best_fit$lambda.M = lamseq[i]
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
        CAMC3_fit(
          y = y,
          X = X,
          J = best_fit$rank.max,
          lambda.M = best_fit$lambda.M,
          learning.rate = hpar$beta$learning.rate,
          beta.iter.max = hpar$beta$prox.iter.max,
          lambda.beta = lambda.beta,
          lambda.a = hpar$laplacian$lambda.a,
          S.a = hpar$laplacian$S.a,
          lambda.b = hpar$laplacian$lambda.b,
          S.b = hpar$laplacian$S.b,
          warm.start = best_fit$fit,
          trace.it = F,
          thresh = thresh,
          maxit = maxit
        )
      
    }
    best_fit$lambda.beta = lambda.beta
    best_fit$lambda.a = hpar$laplacian$lambda.a
    best_fit$lambda.b = hpar$laplacian$lambda.b
    best_fit$learning.rate = hpar$beta$learning.rate
    return(best_fit)
  }
#---
