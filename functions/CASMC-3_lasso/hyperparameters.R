CAMC_Lasso_hparams <-
  list(
    M = list(
      # tuning parameters for lambda
      lambda.factor = 1 / 4,
      lambda.init = NULL,
      n.lambda = 20,
      # tuning parameters for J
      rank.init = 2,
      rank.limit = 30,
      rank.step = 2,
      pct = 0.98,
      early.stopping = 1
    ),
    beta = list(
      # L1 parameters
      learning.rate = "default",
      lambda.max = NULL,
      prox.iter.max = 20,
      n.lambda = 20
    ),
    laplacian = list(
      # laplacian parameters
      lambda.a = 0,
      S.a = NULL,
      lambda.b = 0,
      S.b = NULL
    )
  )
#--------------------------------------------------
# the following function searches for lambda_max for the lasso norm on beta

find_lasso_max_param <- function(dat, hpar=CAMC_Lasso_hparams, verbose=0){
  # 1. Get an upperbound to lambda.beta using the formula:
  
  
  fit0 <- CAMC3_cv_M(
    dat$fit_data$train,
    dat$fit_data$Xq,
    dat$fit_data$valid,
    dat$fit_data$W_valid,
    y = dat$fit_data$Y,
    lambda_beta = 0,
    trace = 0,
    hpar = hpar
  )
  
  E <- dat$fit_data$train - fit0$fit$u %*% (fit0$fit$d * t(fit0$fit$v)) -
    dat$fit_data$Xq %*% fit0$fit$beta
  XtE <- crossprod(dat$fit_data$Xq, E)
  lambda_max =  max(XtE + fit0$fit$beta)
  #------
  # 2. The value above is an upperbound but not the sup. We will look for the sup using
  # line search from 1 -> lambda_max.
  
  # 2.1: Check if the sup is in the first or the second half of the data.
  
  mid_pt <- lambda_max /2
  
  fit1 <- CAMC3_fit(
    dat$fit_data$train,
    dat$fit_data$Xq,
    J = fit0$fit$J,
    lambda.M = fit0$fit$lambda.M,
    lambda.beta = mid_pt,
    trace.it=F
  )
  beta_ratio = sum(fit1$beta ==0) / length(fit1$beta)
  
  if(beta_ratio < 1){
    # 2.2a search in the second half
    old_lambda = lambda_max
    for(lamb in seq(lambda_max, mid_pt, length.out=20)){
      fit2 <- CAMC3_fit(
        dat$fit_data$train,
        dat$fit_data$Xq,
        J = fit0$fit$J,
        lambda.M = fit0$fit$lambda.M,
        lambda.beta = lamb,
        trace.it=F
      )
      beta_ratio = sum(fit2$beta ==0) / length(fit2$beta)
      if (beta_ratio < 1){
        lambda_sup = old_lambda
        break
      }
      old_lambda = lamb
    }
    
  }else{
    # 2.2b search in the first half
    old_lambda = mid_pt
    for(lamb in seq(mid_pt, 0, length.out=20)){
      fit2 <- CAMC3_fit(
        dat$fit_data$train,
        dat$fit_data$Xq,
        J = fit0$fit$J,
        lambda.M = fit0$fit$lambda.M,
        lambda.beta = lamb,
        trace.it=F
      )
      beta_ratio = sum(fit2$beta ==0) / length(fit2$beta)
      if (beta_ratio < 1){
        lambda_sup = old_lambda
        break
      }
      old_lambda = lamb
    }
  }
  if(verbose > 0)
    message(sprintf(
      "The maximum value is %.3f and the sup value is %.3f which is %.1f%% lower.",
      lambda_max, lambda_sup,
      (lambda_sup / lambda_max)*100
    ))
  
  return(lambda_sup)
}