#' Covariate‐Adjusted‐Sparse‐Matrix‐Completion Lasso Cross‐Validation
#'
#' Holdout‐set CV over a grid of \eqn{λ_β}, calling \code{CAMC3_cv_M} under the hood.
#'
#' @param y_train dgCMatrix (Incomplete) of training responses.
#' @param X Covariate matrix.
#' @param y_valid Numeric vector of held‐out responses.
#' @param W_valid Indicator matrix (0 = validation, 1 = training).
#' @param y Optional full dgCMatrix (train + valid) for final refit.
#' @param hpar List of hyperparameters (see \code{CAMC_Lasso_hparams}).
#' @param error_function Function to compute validation error (default = RMSE).
#' @param thresh Convergence tolerance (default = 1e-6).
#' @param maxit Maximum iterations per fit (default = 100).
#' @param trace Integer: 0 = silent, 1 = per‐fit summary, 2 = verbose (passed to inner fits).
#' @param max_cores Maximum parallel workers (default = 8).
#' @param seed Optional RNG seed.
#'
#' @return A list containing the best fit object, its error, hyperparameters, and init params.
#' @export
CAMC_Lasso_cv <- function(
    y_train,
    X,
    y_valid,
    W_valid,
    y               = NULL,
    hpar            = CAMC_Lasso_hparams,
    error_function  = utils$error_metric$rmse,
    thresh          = 1e-6,
    maxit           = 100,
    trace           = 0,
    max_cores       = 8,
    seed            = NULL
) {
  ##-------------------------------------------------------------------------------
  ## 1. Reproducibility & basic checks
  if (!is.null(seed)) set.seed(seed)
  stopifnot(inherits(y_train, "dgCMatrix"))
  
  ##-------------------------------------------------------------------------------
  ## 2. Default β‐learning‐rate if requested
  if (identical(hpar$beta$learning.rate, "default")) {
    hpar$beta$learning.rate <- 1 /
      sqrt(sum((crossprod(X))^2))
  }
  ##-------------------------------------------------------------------------------
  ## 3. Build λ_β grid
  if (is.null(hpar$beta$lambda.max)) {
    nf     <- naive_fit(y_train, X)
    resid  <- y_train - nf$M - X %*% nf$beta
    resid[y_train == 0] <- 0
    hpar$beta$lambda.max <- max(
      (nf$beta / hpar$beta$learning.rate) -
        t(X) %*% resid
    )
  }
  lambda_beta_grid <- seq(
    from       = hpar$beta$lambda.max,
    to         = .Machine$double.eps,
    length.out = hpar$beta$n.lambda
  )
  
  ##-------------------------------------------------------------------------------
  ## 4. Parallel setup
  num_cores <- length(lambda_beta_grid)
  if (num_cores > max_cores) {
    num_cores <- min(max_cores, ceiling(num_cores / 2))
  }
  message("Running on ", num_cores, " cores.")
  
  ##-------------------------------------------------------------------------------
  ## 5. CV loop: fit CAMC3_cv_M for each λ_β in parallel
  inner_trace <- (trace >= 2)
  results     <- parallel::mclapply(
    lambda_beta_grid,
    function(lambda_beta) {
      tryCatch(
        CAMC3_cv_M(
          y_train        = y_train,
          X              = X,
          y_valid        = y_valid,
          W_valid        = W_valid,
          y              = y,
          hpar           = hpar,
          lambda_beta    = lambda_beta,
          error_function = error_function,
          thresh         = thresh,
          maxit          = maxit,
          trace          = inner_trace,
          seed           = seed
        ),
        error = function(e) {
          list(
            error_message = e,
            error         = Inf,
            lambda_beta   = lambda_beta
          )
        }
      )
    },
    mc.cores   = num_cores,
    mc.cleanup = TRUE
  )
  
  ##-------------------------------------------------------------------------------
  ## 6. Report any fit errors
  for (res in results) {
    if (!is.null(res$error_message)) {
      message(
        "Error at λ_β=", round(res$lambda_beta, 3), 
        ": ", res$error_message
      )
    }
  }
  
  ##-------------------------------------------------------------------------------
  ## 7. Select the best fit
  errors    <- vapply(results, `[[`, numeric(1), "error")
  best_idx  <- which.min(errors)
  best_fit  <- results[[best_idx]]
  
  ##-------------------------------------------------------------------------------
  ## 8. Optional printing
  if (trace >= 1) {
    for (res in results) {
      message(sprintf(
        "<< λ_β=%.4g | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
        res$lambda_beta,
        res$error,
        res$fit$n_iter,
        res$fit$J,
        res$fit$lambda.M
      ))
    }
    message(sprintf(
      "<< Best fit >> λ_β=%.4g | err=%.5f | iters=%d | rank_M=%d | λ_M=%.4g >>",
      best_fit$lambda_beta,
      best_fit$error,
      best_fit$fit$n_iter,
      best_fit$fit$J,
      best_fit$fit$lambda.M
    ))
  }
  
  ##-------------------------------------------------------------------------------
  ## 9. Attach metadata & return
  best_fit$hparams   <- data.frame(
    lambda_beta   = best_fit$lambda_beta,
    lambda_M      = best_fit$fit$lambda.M,
    learning_rate = hpar$beta$learning.rate,
    rank_M        = best_fit$fit$J
  )
  best_fit$init_hpar <- hpar
  
  best_fit
}