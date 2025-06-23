Mao_Bixi_Wrapper <-
 function(dat,
          lambda.1_grid = seq(0, 1, length = 20),
          lambda.2_grid = seq(0.9, 0.1, length = 20),
          alpha_grid = c(1),
          ncores = 1,
          # keep it > 1
          n_folds = 5,
          weight_function = Mao_weights$uniform,
          ...) {
  start_time = Sys.time()
  fiti <- Mao.cv(
   Y = dat$Y,
   X = dat$X,
   W = dat$masks$tr_val,
   n_folds = n_folds,
   lambda.1_grid = lambda.1_grid,
   lambda.2_grid = lambda.2_grid,
   alpha_grid = alpha_grid,
   seed = 2023,
   numCores = ncores,
   n1n2_optimized = TRUE,
   test_error = error_metric$rmse,
   theta_estimator = weight_function,
   sequential = FALSE
  )
  
  fit. <- fiti$fit
  
  fit.$Xbeta = dat$X %*% fit.$beta
  
  results = list(model = "Mao")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fiti$best_parameters$lambda.2
  results$lambda.beta = fiti$best_parameters$lambda.1
  
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
    t() |>
    as.data.frame() |>
    mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
      sum(x != 0) / length(x))) |>
    `rownames<-` (colnames(dat$X)) ->
    results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }

SImpute_Bixi_Wrapper <- function(dat, ...) {
 start_time = Sys.time()
 fit. <- simpute.cv(
  Y_train = dat$Y,
  y_valid = dat$splits$valid@x,
  W_valid = dat$masks$valid,
  y = dat$Y,
  n.lambda = 20,
  trace = FALSE,
  print.best = FALSE,
  tol = 5,
  thresh = 1e-6,
  rank.init = 2,
  rank.limit = 30,
  rank.step = 2,
  maxit = 600,
  seed = NULL
 )
 results = list(model = "SoftImpute")
 results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
 results$lambda.beta = NA
 results$lambda.M = fit.$lambda
 results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
 
 results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
 results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
 results$rank_M = qr(fit.$estimates)$rank
 results$rank_beta = NA
 results$sparse_prop = NA
 NA ->
   results$cov_summaries
 
 Y = dat$depart[dat$masks$obs == 1]
 M = fit.$estimates[dat$masks$obs == 1]
 resids = Y - M
 
 TSS = sum((Y - mean(Y)) ^ 2)
 ESS_M = sum((M - mean(Y)) ^ 2)
 RSS <- sum(resids ^ 2)
 
 results$Prop_explained_xbeta = NA
 results$Prop_explained_M = var(M) / var(Y)
 results$Prop_unexplained <-  1 - results$Prop_explained_M 
 
 results
}


SImpute_Bixi_likl <- function(dat, ...) {
  fit. <- simpute.cv(
    Y_train = dat$Y,
    y_valid = dat$splits$valid@x,
    W_valid = dat$masks$valid,
    y = dat$Y,
    n.lambda = 20,
    trace = FALSE,
    print.best = FALSE,
    tol = 5,
    thresh = 1e-6,
    rank.init = 2,
    rank.limit = 30,
    rank.step = 2,
    maxit = 600,
    seed = NULL
  )
  
  
  Y = dat$depart[dat$masks$test == 0]
  M = fit.$estimates[dat$masks$test == 0]
  resids = Y - M
  sigmasq = var(resids)
  log.likli = - length(resids)/2 * log(2 * pi * sigmasq) - 1 / (2*sigmasq) * sum(resids^2)
  return(log.likli)
  
}


#----------------------------------------------------------------------------------
CAMC_1_Bixi_Wrapper <-
 function(dat,
          max_cores = 20,
          maxit = 300,
          ...) {
  start_time = Sys.time()
  
  fiti <- CAMC1_cv(
   y_train = dat$splits$train,
   X = dat$X,
   y_valid = dat$splits$valid,
   W_valid = dat$masks$valid,
   y = dat$depart,
   error_function = error_metric$rmse,
   lambda.factor = 1 / 4,
   lambda.init = NULL,
   n.lambda = 20,
   rank.init = 2,
   rank.limit = 30,
   rank.step = 2,
   pct = 0.98,
   lambda.a = 0,
   S.a = NULL,
   lambda.b = 0,
   S.b = NULL,
   early.stopping = 1,
   thresh = 1e-6,
   maxit = maxit,
   trace = F,
   print.best = TRUE,
   quiet = FALSE,
   warm = NULL,
   r_min = 0,
   track = F,
   max_cores = max_cores,
   seed = NULL
  )
  
  fit. = fiti$fit
  # get estimates and validate
  fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
  fit.$Xbeta = dat$X %*% fit.$beta
  fit.$estimates = fit.$M + fit.$Xbeta
  
  results = list(model = "CAMC-1")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fit.$lambda
  results$lambda.beta = NA
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
    t() |>
    as.data.frame() |>
    mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
      sum(x != 0) / length(x))) |>
    `rownames<-` (colnames(dat$X)) ->
    results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }

#--------------------------------------------------------------------------------------
CAMC_0_Bixi_Wrapper <-
 function(dat,
          max_cores = 20,
          maxit = 300,
          ...) {
  start_time = Sys.time()
  
  fiti <- CAMC0_cv(
   y_train = dat$splits$train,
   X = dat$X,
   y_valid = dat$splits$valid@x,
   W_valid = dat$masks$valid,
   #y = dat$depart,
   error_function = error_metric$rmse,
   lambda.factor = 1 / 4,
   lambda.init = NULL,
   n.lambda = 20,
   rank.init = 2,
   rank.limit = 30,
   rank.step = 2,
   pct = 0.98,
   lambda.a = 0,
   S.a = NULL,
   lambda.b = 0,
   S.b = NULL,
   early.stopping = 1,
   thresh = 1e-6,
   maxit = maxit,
   trace = FALSE,
   print.best = F,
   quiet = FALSE,
   warm = NULL,
   lambda.beta.grid = "default",
   track = F,
   max_cores = max_cores,
   seed = NULL
  )
  
  fit. = fiti$fit
  # get estimates and validate
  fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
  fit.$Xbeta = dat$X %*% fit.$beta
  fit.$estimates = fit.$M + fit.$Xbeta
  
  results = list(model = "CAMC-0")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fit.$lambda
  results$lambda.beta = fiti$lambda.beta
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
    t() |>
    as.data.frame() |>
    mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
      sum(x != 0) / length(x))) |>
    `rownames<-` (colnames(dat$X)) ->
    results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }
#-------


CAMC_2_Bixi_Wrapper <-
 function(dat,
          max_cores = 20,
          maxit = 300,
          ...) {
  start_time = Sys.time()
  
  fiti <- CAMC2_cv2(
   y_train = dat$splits$train,
   X = dat$X,
   y_valid = dat$splits$valid@x,
   W_valid = dat$masks$valid,
   #y = dat$depart,
   error_function = error_metric$rmse,
   warm = NULL,
   M_cv_param = list(
     rank.init = 2,
     rank.limit = 30,
     rank.step = 2,
     pct = 0.98,
     lambda.factor = 1/4,
     lambda.init = NULL,
     n.lambda = 20, 
     early.stopping = 1
   ),
   beta_cv_param = list(
     rank.init = 2,
     rank.limit = qr(dat$X)$rank,
     rank.step = 2,
     pct = 0.98,
     lambda.multi.factor = 20,
     lambda.init = NULL,
     n.lambda = 20, 
     early.stopping = 1
   ),
   quiet = T,
   trace = F,
   track = F,
   step3 = T,
   use_warmstart = TRUE,
   seed = NULL,
  )
  
  fit. = fiti$fit
  # get estimates and validate
  fit.$M = unsvd(fit.) 
  fit.$beta = unsvd(fit.$beta)
  fit.$Xbeta = dat$X %*% fit.$beta
  fit.$estimates = fit.$M + fit.$Xbeta
  
  results = list(model = "CAMC-2")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fiti$hparams$lambda.M
  results$lambda.beta = fiti$hparams$lambda.beta
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
    t() |>
    as.data.frame() |>
    mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
      sum(x != 0) / length(x))) |>
    `rownames<-` (colnames(dat$X)) ->
    results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }
#----------------------------------------------------
CAMC_3a_Bixi_Wrapper <-
 function(dat,
          max_cores = 20,
          maxit = 300,
          ...) {
  start_time = Sys.time()
  learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X) ^ 2))
  
  
  fiti <- CAMC3_cv_beta(
   y_train = dat$splits$train,
   X = dat$X,
   y_valid = dat$splits$valid@x,
   W_valid = dat$masks$valid,
   # y = to_incomplete(dat$Y),
   trace = 0,
   print.best = T,
   warm = NULL,
   quiet = F,
   learning.rate = learning_rate,
   early.stopping = 1,
   lambda.beta.grid = seq(0, 10, length.out = 20),
   max_cores = max_cores
  )
  
  fit. = fiti$fit
  # get estimates and validate
  fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
  fit.$Xbeta = dat$X %*% fit.$beta
  fit.$estimates = fit.$M + fit.$Xbeta
  
  results = list(model = "CAMC-3a")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fiti$hparams$lambda.M
  results$lambda.beta = fiti$hparams$lambda.beta
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
   t() |>
   as.data.frame() |>
   mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
    sum(x != 0) / length(x))) |>
   `rownames<-` (colnames(dat$X)) ->
   results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }
#------
CAMC_3b_Bixi_Wrapper <-
 function(dat,
          max_cores = 20,
          maxit = 300,
          ...) {
  start_time = Sys.time()
  learning_rate = 1 / sqrt(sum((t(dat$X) %*% dat$X) ^ 2))
  fiti <- CAMC3_kfold(
   Y = dat$Y,
   X = dat$X,
   obs_mask = dat$masks$tr_val,
   n_folds = 5,
   trace = 2,
   print.best = T,
   warm = NULL,
   beta.iter.max = 30,
   quiet = F,
   learning.rate = learning_rate,
   early.stopping = 1,
   lambda.beta.grid = seq(0, 2.5, length.out = 10),
   max_cores = max_cores
  )
  
  fit. = fiti$fit
  # get estimates and validate 
  fit.$M = fit.$u %*% (fit.$d * t(fit.$v))
  fit.$Xbeta = dat$X %*% fit.$beta 
  fit.$estimates = fit.$M + fit.$Xbeta
  
  results = list(model = "CAMC-3b")
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  results$lambda.M = fiti$hparams$lambda.M
  results$lambda.beta = fiti$hparams$lambda.beta
  results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
  results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
  results$rank_M = sum(fit.$d > 0)
  results$rank_beta = qr(fit.$beta)$rank
  results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
  apply(fit.$beta, 1, summary) |> as.data.frame() |>
    t() |>
    as.data.frame() |> 
    mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
      sum(x != 0) / length(x))) |>
    `rownames<-` (colnames(dat$X)) ->
    results$cov_summaries
  
  Y = dat$depart[dat$masks$obs == 1]
  xbeta = fit.$Xbeta[dat$masks$obs == 1]
  M = fit.$M[dat$masks$obs == 1]
  resids = Y - xbeta - M
  
  TSS = sum((Y - mean(Y)) ^ 2)
  ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
  ESS_M = sum((M - mean(Y)) ^ 2)
  RSS <- sum(resids ^ 2)
  
  results$Prop_explained_xbeta = var(xbeta) / var(Y)
  results$Prop_explained_M = var(M) / var(Y)
  results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
  
  Y = dat$depart[dat$masks$test == 0]
  xbeta = fit.$Xbeta[dat$masks$test == 0]
  M = fit.$M[dat$masks$test == 0]
  resids = Y - xbeta - M
  
  
  sigmasq = var(resids)
  log.likli = (- length(resids)/2) * log(2 * pi * sigmasq) - (1 / (2*sigmasq)) * sum(resids^2)
  Rsq.mcf = 1 - log.likli / SImpute_Bixi_likl(dat)
  results$Prop_explained_xbeta <- Rsq.mcf
  
  results
 }
#------





#----------------------------------------------
Naive_Bixi_Wrapper <- function(dat, ...) {
 start_time = Sys.time()
 fiti <- naive_fit(dat$Y, dat$X)
 results = list(model = "Naive")
 fit. <- fiti
 fit.$Xbeta = dat$X %*% fit.$beta
 
 results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
 results$lambda.M = NA
 results$lambda.beta = NA
 results$error.test = test_error(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
 results$corr.test = cor(fit.$estimates[dat$masks$test == 0], dat$splits$test@x)
 results$error.train = test_error(fit.$estimates[dat$masks$tr_val != 0], dat$depart@x)
 results$rank_M = sum(fit.$d > 0)
 results$rank_beta = qr(fit.$beta)$rank
 results$sparse_prop = sum(fit.$beta == 0) / length(fit.$beta)
 apply(fit.$beta, 1, summary) |> as.data.frame() |>
   t() |>
   as.data.frame() |>
   mutate(prop_non_zero = apply(fit.$beta, 1, function(x)
     sum(x != 0) / length(x))) |>
   `rownames<-` (colnames(dat$X)) ->
   results$cov_summaries
 
 Y = dat$depart[dat$masks$obs == 1]
 xbeta = fit.$Xbeta[dat$masks$obs == 1]
 M = fit.$M[dat$masks$obs == 1]
 resids = Y - xbeta - M
 
 TSS = sum((Y - mean(Y)) ^ 2)
 ESS_xbeta = sum((xbeta - mean(Y)) ^ 2)
 ESS_M = sum((M - mean(Y)) ^ 2)
 RSS <- sum(resids ^ 2)
 
 results$Prop_explained_xbeta = var(xbeta) / var(Y)
 results$Prop_explained_M = var(M) / var(Y)
 results$Prop_unexplained <-  1 - results$Prop_explained_M -  results$Prop_explained_xbeta 
 
 
 
 results
}
