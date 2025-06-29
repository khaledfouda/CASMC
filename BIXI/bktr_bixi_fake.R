setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")

model.dat <- load_bixi_dat(transpose = T, scale_response = F)$model
#------------------------------------------------------------------
print_performance <-
  function(inp,
           outp,
           test_error = error_metric$rmse,
           mao = FALSE,
           name = "",
           showres = T,
           rdig = 4) {
    error_function_name <-
      names(which(sapply(error_metric, identical, test_error)))
    
    if (mao) {
      M = outp$M
      O = outp$estimates
    } else{
      M = outp$u %*% (outp$d * t(outp$v))
      O = M + outp$X %*% outp$beta
    }
    
    error_xbeta <- error_M <- rank_beta <- rank_M <- NA
    if (!is.null(outp$beta))
      error_xbeta <-
      tryCatch(round(test_error(outp$X %*% outp$beta, inp$X %*% inp$beta), rdig),error=function(e)NA)
    if (!is.null(M))
      error_M <- 
      tryCatch(round(test_error(M, inp$M), rdig),error=function(e)NA)
    error_test <-
      tryCatch(round(test_error(O[inp$W == 0], inp$O[inp$W == 0]), rdig),error=function(e)NA)
    error_train <-
      tryCatch(round(test_error(O[inp$W != 0], inp$O[inp$W != 0]), rdig),error=function(e)NA)
    error_Y <-
      tryCatch(round(test_error(O[inp$W != 0], inp$Y[inp$W != 0]), rdig),error=function(e)NA)
    
    if (!is.null(outp$beta))
      rank_beta <- qr(outp$beta)$rank
    if (!is.null(M))
      rank_M <- qr(M)$rank
    
    result_df <- data.frame(
      Metric = c(paste0(error_function_name, "(rank)")),
      model = name,
      XBeta = tryCatch(sprintf(paste0("%.", rdig, "f(%2d)"), error_xbeta, rank_beta),error=function(e)NA),
      M = tryCatch(sprintf(paste0("%.", rdig, "f(%2d)"), error_M, rank_M),error=function(e)NA),
      test = tryCatch(sprintf(paste0("%.", rdig, "f"), error_test),error=function(e)NA),
      train = tryCatch(sprintf(paste0("%.", rdig, "f"), error_train),error=function(e)NA),
      Y = tryCatch(sprintf(paste0("%.", rdig, "f"), error_Y),error=function(e)NA)
    )
    
    if (showres) {
      print(knitr::kable(result_df, format = "simple"))
    } else
      return(result_df)
  }







#-------------------------------------------------------
generate_fake_bixi <- function(dat,
                               r = 5,
                               missing_prob = 0.9,
                               informative_cov_prop = 1,
                               prepare_for_fitting = TRUE,
                               beta_is_cor = FALSE,
                               seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  X <- dat$X
  Y <- naive_MC(as.matrix(dat$depart))
  k = ncol(X)
  n = nrow(X)
  m = ncol(dat$depart)
  if (beta_is_cor) {
    beta <- matrix(0, k, m)
    for (i in 1:k) {
      for (j in 1:m)
        beta[i, j] <- cor(X[, i], Y[, j])
    }
    
  } else{
    beta_means <- runif(k, 1, 3) * sample(c(-1, 1), k, TRUE)
    beta_vars <- runif(k, 0.5, 1) ^ 2
    beta <-
      t(mvrnorm(m, beta_means, diag(beta_vars, k, k)))
  }
  U <- matrix(runif(n * r), ncol = r)
  V <- matrix(runif(r * m), nrow = r)
  P_X = X %*% solve(t(X) %*% X) %*% t(X)
  P_bar_X = diag(1, n, n) - P_X
  
  M <- P_bar_X %*% U %*% V
  
  W <- matrix(rbinom(n * m, 1, (1 - missing_prob)) , nrow = n)
  ncov_to_keep = round(informative_cov_prop * k)
  
  if (ncov_to_keep <= 0) {
    beta <- matrix(0, k, ncol = m)
  } else if (ncov_to_keep < k) {
    sampled_covars_to_remove <-
      sample(1:k, k - ncov_to_keep, replace = FALSE)
    beta[sampled_covars_to_remove, ] <- 0
  }
  O <- X %*% beta +  M
  E <-
    matrix(rnorm(n * m, mean = 0, sd = 1), ncol = m)
  Y <- (O + E) * W
  rank <- qr(O)$rank
  
  fit_data <- NULL
  if (prepare_for_fitting) {
    fit_data <- list()
    fit_data$W_valid <- matrix.split.train.test(W, testp = 0.2)
    train = (Y * fit_data$W_valid)
    fit_data$train = to_incomplete(train)
    fit_data$valid = Y[fit_data$W_valid == 0]
    fit_data$Y = to_incomplete(Y)
  }
  
  return(list(
    O = O,
    W = W,
    X = X,
    Y = Y,
    beta = beta,
    M = M,
    fit_data = fit_data,
    rank = rank
  ))
}

#----------------------------------------------------------------------------------


bixi.dat <-
  load_bixi_dat(transpose = F,
                scale_covariates = F,
                seed = 2023)$model
model.dat <-
  generate_fake_bixi(
    bixi.dat,
    r =  15,
    missing_prob =  0.5,
    beta_is_cor = TRUE,
    informative_cov_prop = 0.7,
    seed = 2023
  )
print(paste("True rank of beta is", qr(model.dat$beta)$rank))
print(paste("True rank of M is", qr(model.dat$M)$rank ))

all_res = data.frame()
case = 1
for (case in 0:4) {
  if (case == 0) {
    X <- model.dat$X
  } else if (case == 1) {
    X <- model.dat$X |>
      scalers("minmax")
  } else if (case == 2) {
    X <- model.dat$X |>
      scalers("minmax") |>
      remove_collinear_cols(thresh = 0.7)
  } else if (case == 3) {
    X <- remove_collinear_cols(model.dat$X, 0.7)
  } else if (case == 4) {
    X <- reduced_hat_decomp(model.dat$X, .01, 0.98)$X
  }
  
  results <- CAMC_var_selection(
    y_train = model.dat$fit_data$train,
    y_valid = model.dat$fit_data$valid,
    Y = model.dat$Y,
    X = X,
    W_valid = model.dat$fit_data$W_valid,
    track = F,
    return_best_model = TRUE,
  )
  
  results$res |>
    select(-variables) |>
    arrange(validation_error) |>
    kable() |>
    print()
  fiti <- results$fit$fit
  fiti$X <- results$fit$X
  
  all_res <- rbind(all_res,
                   print_performance(
                     model.dat,
                     fiti,
                     error_metric$rmse,
                     F,
                     paste0("CAMC(", case, ")"),
                     F,
                     3
                   ))
  print(all_res)
}


simpute_fit <- simpute.cv(
  Y_train = as.matrix(model.dat$fit_data$train),
  y_valid = model.dat$fit_data$valid,
  W_valid = model.dat$fit_data$W_valid,
  y = model.dat$Y,
  n.lambda = 20,
  trace = FALSE,
  print.best = TRUE,
  tol = 5,
  thresh = 1e-6,
  rank.init = 2,
  rank.limit = 30,
  rank.step = 2,
  maxit = 300,
  seed= 2023
  
)


all_res  |>
  rbind(print_performance(model.dat, simpute_fit, error_metric$rmse, TRUE, "SoftImpute", F, 3)) |> 
  arrange(XBeta) |>
  mutate(XBeta = paste0("[", 1:nrow(all_res), "]", XBeta)) |>
  arrange(M) |>
  mutate(M = paste0("[", 1:nrow(all_res), "]", M)) |>
  arrange(Y) |>
  mutate(Y = paste0("[", 1:nrow(all_res), "]", Y)) |>
  dplyr::select(-Y) |>
  arrange(train) |>
  mutate(train = paste0("[", 1:nrow(all_res), "]", train)) |>
  arrange(test) |>
  mutate(test = paste0("[", 1:nrow(all_res), "]", test))  ->
  all_res


kable(all_res, format = "simple") |>  print()


#
# fiti$beta
# apply(fiti$fit$beta, 1, summary) |> print()
# test_error <- error_metric$rmse
# fit1 = best_fit$fit
# sout = best_fit
# # get estimates and validate
# sout$M = unsvd(fit1)
# sout$beta =  fit1$beta
# sout$estimates = sout$M + X %*% (sout$beta)
#
#
# which(model.dat$beta[,1]==0)
#
# fiti$beta[,1]
# model.dat$beta[,1]
#
#
#
#
#
# #----
# # apply models:
# model.dat$X[1,]
# # CAMC
# aresults <- list()
# i = 1
#
#
# sparm = list(NULL,NULL,"")
# for(sparm in list(list(NULL,NULL,""),
#                   list(model.dat$similarity.A, NULL, "_with_time_laplacian")
#                   #list(NULL, model.dat$sim_col),
#                   #list(model.dat$sim_row, model.dat$sim_col)
#                   )
#     ){
#
#   #X = cbind(X[,1:4]^2)
#   #------------------------
#   # X <- cbind(
#   #   #X,
#   #   scale(X[,1:3])^2,
#   #   log(X[,1, drop=FALSE]+abs(min(X[,1]))+1)
#   #   #matrix(X[,4]>0,ncol= 1)
#   # )
#
#   #---------------------
#   corr <- cor(model.dat$X) |> round(2); corr[corr<0.5] = diag(corr) = corr[upper.tri(corr)] = NA; corr
#   #------------------------------------------------------------------------
#   for(b in 1:5){
#   #model.dat <- load_bixi_dat(transpose = T, scale_response = T, seed=b)$model
#   start_time = Sys.time()
#
#   bixi.dat <- load_bixi_dat(transpose = F, scale_covariates = F)$model
#   #bixi.dat$X <- apply(bixi.dat$X, 2, function(x) (x - min(x)) / (max(x) - min(x)))
#   model.dat <- generate_fake_bixi(bixi.dat, 15, 0.3,seed = 2023)
#   model.dat$Xo <- model.dat$X
#   model.dat$X <- apply(model.dat$X, 2, function(x) (x - min(x)) / (max(x) - min(x)))
#   # model.dat <- generate_simulation_rows(300,460,5,7,0.9,FALSE,FALSE,.7,T)
#   bixi.dat$X[1:3,] |> t() |> as.data.frame() |>  mutate(id=1:ncol(bixi.dat$X)) |>  kable()
#   X <-  model.dat$Xo[,-c(1,5,9,12,16)]
#   X <- X[,-c(8,11)]
#   X_r <- reduced_hat_decomp(X, .02, 0.99)
#   X_r$rank
#   X <- X_r$X
#   X[1,]
#   best_fit = CAMC_cv_rank(
#     y_train = model.dat$fit_data$train,
#     X = X,#[,1:5, drop=FALSE],
#
#     y_valid = model.dat$fit_data$valid,
#     W_valid = model.dat$fit_data$W_valid,
#     y = model.dat$fit_data$Y,
#     trace = F,
#     max_cores = 100,
#     thresh = 1e-6,
#     lambda.a = 0.01,
#     S.a = sparm[[1]],
#     lambda.b = 0.2,
#     S.b = sparm[[2]],
#     #n.lambda = n.lambda,
#     #rank.limit = 30,
#     maxit = 600,
#     r_min = 0,
#     rank.init = 10,
#     rank.step = 4,
#     print.best = TRUE,
#     seed = 2023,
#     track  = T
#   )
#   best_fit$rank_M
#   apply(best_fit$fit$beta, 1, summary) |> print()
#   test_error <- error_metric$rmse
#   fit1 = best_fit$fit
#   sout = best_fit
#   # get estimates and validate
#   sout$M = unsvd(fit1)
#   sout$beta =  fit1$beta
#   sout$estimates = sout$M + X %*% (sout$beta)
#
#   hist(sout$beta[1,])
#
#   plot(1:length(sout$estimates), as.vector(sout$estimates))
#
#   #plot(1:ncol(sout$beta), sout$beta[3,])
#
#   if(b > 1) old_results = results
#   results = list(model = paste0("CAMC_rank_all",sparm[[3]]))
#   results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#   #results$lambda.1 = NA#sout$lambda.beta |> round(3)
#   #results$lambda.2 = sout$lambda |> round(3)
#   results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
#                                   model.dat$splits$test@x) |> round(5)
#   results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                    model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#     round(5)
#   results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
#                                    model.dat$splits$valid@x) |> round(5)
#
#   results$rank = qr(sout$estimates)$rank
#   if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
#   }
#   results[-1] = mapply(function(x)x/5, results[-1])
#
#
#   aresults[[i]] <- results
#   i = i +1
#   #------------------------------
#
#
# X <- X[,1:3]
#   for(b in 1:5){
#     model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
#
# best_fit = CAMC_cv_rank(
#  y_train = model.dat$splits$train,
#  X = X,
#  y_valid = model.dat$splits$valid@x,
#  W_valid = model.dat$masks$valid ,
#  y = model.dat$depart,
#  trace = F,
#  max_cores = 30,
#  thresh = 1e-6,
#  lambda.a = 0.01,
#  S.a = sparm[[1]],
#  lambda.b = 0.2,
#  S.b = sparm[[2]],
#  #n.lambda = n.lambda,
#  #rank.limit = rank.limit,
#  maxit = 200,
#  #rank.step = rank.step,
#  print.best = TRUE,
#  seed = 2023,
#  track  = T
# )
# test_error <- error_metric$rmse
# fit1 = best_fit$fit
# sout = best_fit
# # get estimates and validate
# sout$M = unsvd(fit1)
# sout$beta =  fit1$beta
# apply(sout$beta, 1, summary) |> print()
# sout$estimates = sout$M + X %*% (sout$beta)
# if(b > 1) old_results = results
#
# results = list(model = paste0("CAMC_rank_1:3_",sparm[[3]]))
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
# #results$lambda.1 = NA#sout$lambda.beta |> round(3)
# #results$lambda.2 = sout$lambda |> round(3)
# results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
#                                 model.dat$splits$test@x) |> round(5)
# results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                 model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#  round(5)
# results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
#                                 model.dat$splits$valid@x) |> round(5)
#
# results$rank = qr(sout$estimates)$rank
# if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
#   }
#   results[-1] = mapply(function(x)x/5, results[-1])
#
#
# aresults[[i]] <- results
# i = i +1
# #-----------------------------
#
# X <- X[,2, drop=FALSE]
#
# for(b in 1:5){
#   model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
#
# best_fit = CAMC_cv_rank(
#   y_train = model.dat$splits$train,
#   X = X,
#   y_valid = model.dat$splits$valid@x,
#   W_valid = model.dat$masks$valid ,
#   y = model.dat$depart,
#   trace = F,
#   max_cores = 30,
#   thresh = 1e-6,
#   lambda.a = 0.01,
#   S.a = sparm[[1]],
#   lambda.b = 0.2,
#   S.b = sparm[[2]],
#   #n.lambda = n.lambda,
#   #rank.limit = rank.limit,
#   maxit = 200,
#   #rank.step = rank.step,
#   print.best = TRUE,
#   seed = 2023,
#   track  = T
# )
# test_error <- error_metric$rmse
# fit1 = best_fit$fit
# sout = best_fit
# # get estimates and validate
# sout$M = unsvd(fit1)
# sout$beta =  fit1$beta
# sout$estimates = sout$M + X %*% (sout$beta)
# if(b > 1) old_results = results
#
# results = list(model = paste0("CAMC_rank_2_",sparm[[3]]))
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
# #results$lambda.1 = NA#sout$lambda.beta |> round(3)
# #results$lambda.2 = sout$lambda |> round(3)
# results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
#                                 model.dat$splits$test@x) |> round(5)
# results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                  model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#   round(5)
# results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
#                                  model.dat$splits$valid@x) |> round(5)
#
# results$rank = qr(sout$estimates)$rank
# if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
# }
# results[-1] = mapply(function(x)x/5, results[-1])
#
#
# aresults[[i]] <- results
# i = i +1
# #------------------------------
# X <- model.dat$X |> scale()
# X = cbind(X[,1:4]^2)
#
# for(b in 1:5){
#   model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
#
# best_fit = CAMC_cv_L2(
#   y_train = model.dat$splits$train,
#   X = X,
#   y_valid = model.dat$splits$valid@x,
#   W_valid = model.dat$masks$valid ,
#   y = model.dat$depart,
#   trace = F,
#   max_cores = 30,
#   thresh = 1e-6,
#   lambda.a = 0.01,
#   S.a = sparm[[1]],
#   lambda.b = 0.2,
#   S.b = sparm[[2]],
#   #n.lambda = n.lambda,
#   #rank.limit = rank.limit,
#   maxit = 200,
#   #rank.step = rank.step,
#   print.best = TRUE,
#   seed = 2023,
#   track = T
# )
# test_error <- error_metric$rmse
# fit1 = best_fit$fit
# sout = best_fit
# # get estimates and validate
# sout$M = unsvd(fit1)
# sout$beta =  fit1$beta
# apply(sout$beta, 1, summary) |> print()
# sout$estimates = sout$M + X %*% (sout$beta)
# if(b > 1) old_results = results
#
# results = list(model = paste0("CAMC_rank_L2_all_",sparm[[3]]))
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
# #results$lambda.1 = NA#sout$lambda.beta |> round(3)
# #results$lambda.2 = sout$lambda |> round(3)
# results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
#                                 model.dat$splits$test@x) |> round(5)
# results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                  model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#   round(5)
# results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
#                                  model.dat$splits$valid@x) |> round(5)
#
# results$rank = qr(sout$estimates)$rank
# if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
# }
# results[-1] = mapply(function(x)x/5, results[-1])
#
#
# aresults[[i]] <- results
# i = i +1
#
# }
# #-------------------------------------------------------------------------------
# # Soft Impute
# for(b in 1:5){
#   model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
# sout <- simpute.cv(
#  Y_train = as.matrix(model.dat$splits$train),
#  y_valid = model.dat$splits$valid@x,
#  W_valid = model.dat$masks$valid,
#  y = as.matrix(model.dat$depart),
#  trace = T,
#  rank.limit = 30,
#  tol = 10,
#  print.best = FALSE,
#  rank.step = 2,
#  maxit = 700
# )
#
# if(b > 1) old_results = results
#
# results = list(model = "SoftImpute")
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
# #results$lambda.1 = NA
# #results$lambda.2 = sout$lambda |> round(3)
# results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
#                                 model.dat$splits$test@x) |> round(5)
# results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                  model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#  round(5)
# results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
#                                  model.dat$splits$valid@x) |> round(5)
#
# results$rank = sout$rank_M
# if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
# }
# results[-1] = mapply(function(x)x/5, results[-1])
#
# aresults[[i]] <- results
# i = i +1
# #-----------------------------------------------------------------------------------
# for(b in 1:5){
#   model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
# estimates = naive_MC(as.matrix(model.dat$splits$train))
# if(b > 1) old_results = results
#
# results = list(model = "Naive")
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
# #results$lambda.1 = NA
# #results$lambda.2 = NA
# results$error.test = test_error(estimates[model.dat$masks$test == 0],
#                                 model.dat$splits$test@x) |> round(5)
# results$error.train = test_error(estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
#                                  model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |>
#  round(5)
# results$error.valid = test_error(estimates[model.dat$masks$valid == 0],
#                                  model.dat$splits$valid@x) |> round(5)
# results$rank = qr(estimates)$rank
# if(b > 1)
#   results[-1] = mapply(sum, old_results[-1], results[-1])
# }
# results[-1] = mapply(function(x)x/5, results[-1])
#
# aresults[[i]] <- results
# i = i +1
# #---------------------------------------------------------------------------------
# do.call(rbind, lapply(aresults, function(x) data.frame(t(unlist(x))))) |> kable()
#
