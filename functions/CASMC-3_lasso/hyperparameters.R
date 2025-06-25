#' Default hyperparameters for CAMC Lasso cross‐validation
#'
#' A nested list of tuning parameters used by \code{CAMC_Lasso_cv} and \code{CAMC3_cv_M}.
#'
#' @format A list with components:
#' \describe{
#'   \item{M}{Parameters for the low‐rank penalty and rank control}
#'   \item{beta}{Parameters for the L1 (lasso) penalty and β updates}
#'   \item{laplacian}{Parameters for optional Laplacian regularization}
#' }
#' @export
CAMC_Lasso_hparams <- list(
  M = list(
    lambda_init    = NULL, # = lambda_{M,min}
    lambda_factor  = 1 / 4, # ignored if lambda_init is provided
    n.lambda       = 20, # seq(lamba.init, 0, length=n.lambda)
    rank.init      = 2, 
    rank.limit     = 30,
    rank.step      = 2,
    early.stopping = 1
  ),
  beta = list(
    lambda_max     = NULL, # if NULL then it's computed. use NULL.
    n.lambda       = 20
  ),
  laplacian = list(
    lambda_a       = 0,
    S.a            = NULL,
    lambda_b       = 0,
    S.b            = NULL
  )
)

#' Find supremum λ for the Lasso penalty on β
#'
#' Computes the smallest \eqn{λ_β} such that all β coefficients are zero,
#' by first obtaining an upper bound and then performing a line search.
#'
#' @param data A list containing \code{fit_data} with elements:
#'   \code{train}, \code{Xq}, \code{valid}, \code{W_valid}, and \code{Y}.
#' @param hparams Hyperparameter list, typically \code{CAMC_Lasso_hparams}.
#' @param verbose Integer; if > 0, prints diagnostic messages (default = 0).
#' @return Numeric scalar: the supremum value of \eqn{λ_β}.
#' @export
find_lasso_max_param <- function(
    y_train,
    X,
    y_valid = NULL,
    W_valid = NULL,
    y       = NULL,
    hparams = CAMC_Lasso_hparams,
    maxit = 100,
    verbose = 0
) {
  ## 1. Initial fit at λ_β = 0
  if(is.null(y_valid) | is.null(W_valid)){
    fit0 <- list()
    fit0$fit <- CAMC3_fit(
      y            = y_train,
      X            = X,
      J            = 10,
      lambda_M     = 1,
      lambda_beta  = 0,
      maxit        = maxit, 
      trace        = FALSE
    )
  }else{
    fit0 <- CAMC3_cv_M(
      y_train        = y_train,
      X              = X,
      y_valid        = y_valid,
      W_valid        = W_valid,
      y              = y,
      lambda_beta    = 0,
      trace          = 0,
      maxit          = maxit, 
      hpar           = hparams
    )
  }
  
  ## 2. Compute λ_max = max( XᵀE + β )
  residuals   <- y_train -
    fit0$fit$u %*% (fit0$fit$d * t(fit0$fit$v)) -
    X %*% fit0$fit$beta
  xt_resid    <- crossprod(X, residuals)
  lambda_max  <- max(xt_resid + fit0$fit$beta)
  
  ## 3. Line‐search for supremum λ_β
  mid_pt      <- lambda_max / 2
  
  fit1 <- CAMC3_fit(
    y            = y_train,
    X            = X,
    J            = fit0$fit$J,
    lambda_M     = fit0$fit$lambda_M,
    lambda_beta  = mid_pt,
    maxit        = maxit,
    trace        = FALSE
  )
  zero_ratio <- sum(fit1$beta == 0) / length(fit1$beta)
  
  if (zero_ratio < 1) {
    # search in [mid_pt, lambda_max]
    old_lambda <- lambda_max
    for (lam in seq(lambda_max, mid_pt, length.out = 20)) {
      fit2 <- CAMC3_fit(
        y            = y_train,
        X            = X,
        J            = fit0$fit$J,
        lambda_M     = fit0$fit$lambda_M,
        lambda_beta  = lam,
        maxit        = maxit,
        trace        = FALSE
      )
      zero_ratio <- sum(fit2$beta == 0) / length(fit2$beta)
      if (zero_ratio < 1) {
        lambda_sup <- old_lambda
        break
      }
      old_lambda <- lam
    }
  } else {
    # search in [0, mid_pt]
    old_lambda <- mid_pt
    for (lam in seq(mid_pt, 0, length.out = 20)) {
      fit2 <- CAMC3_fit(
        y            = y_train,
        X            = X,
        J            = fit0$fit$J,
        lambda_M     = fit0$fit$lambda_M,
        lambda_beta  = lam,
        maxit        = maxit,
        trace        = FALSE
      )
      zero_ratio <- sum(fit2$beta == 0) / length(fit2$beta)
      if (zero_ratio < 1) {
        lambda_sup <- old_lambda
        break
      }
      old_lambda <- lam
    }
  }
  
  if (verbose > 0) {
    message(sprintf(
      "λ_max = %.3f; λ_sup = %.3f (%.1f%% of λ_max)",
      lambda_max,
      lambda_sup,
      100 * lambda_sup / lambda_max
    ))
  }
  
  lambda_sup
}