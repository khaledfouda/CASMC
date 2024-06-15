setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")

#----------------------------------------------------------------------------------
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
      tryCatch(round(test_error(O[inp$masks$test == 0], inp$splits$test[inp$masks$test==0]), rdig),error=function(e)NA)
    error_train <-
      tryCatch(round(test_error(O[inp$masks$obs != 0], inp$depart[inp$masks$obs != 0]), rdig),error=function(e)NA)
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


#----------------------------------------------
model.dat <-
  load_bixi_dat(transpose = F, scale_response = T, scale_covariates = F,
                testp = 0.3, validp = 0.1, seed=2023)$model

model.dat$X <- model.dat$spatial_simple
all_res = data.frame()
case =1
#for (case in 0:4) {
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
    X <- model.dat$X
    X[,colnames(X)[colnames(X) !="num_metro_stations"]] <- 
      X[,colnames(X)[colnames(X) !="num_metro_stations"]] |> 
    scalers("median") #|>
      #remove_collinear_cols(thresh = 0.7)
  }
  
  CASMC2_cv_beta(
    y_train = model.dat$splits$train,
    X = X,
    y_valid = model.dat$splits$valid@x,
    W_valid = model.dat$masks$valid,
    y = model.dat$depart,
    trace = T,
    quiet = F,
    seed = 2023,
    rank.beta.init = 1,
    early.stopping = 1,
    lambda.beta.grid = "default1"
  ) -> fit_bixi
  
  fit_bixi$hparams
  fit_bixi$fit -> fiti3
  fiti3$beta = as.matrix(fiti3$ub %*% (fiti3$db^2) %*% t(fiti3$vb))
  fiti3$X = X
  fiti3$beta[,1:5]
  print_performance(model.dat, fiti3, error_metric$rmse, F, "SoftImpute", F, 3)
  apply(fiti3$beta, 1, summary) |> as.data.frame() |> round(2) |>
    t() |> as.data.frame() |>  arrange(desc(Median)) |>  kable()
  
  
  length(model.dat$depart@x)
  length(model.dat$splits$train@x)
  length(model.dat$splits$test@x)
  sum(model.dat$masks$valid==0)
  sum(model.dat$masks$test==0)
  #---------------------------------------------------------------------
  
  CASMC3_cv_beta(
    y_train = model.dat$splits$train,
    X = X,
    y_valid = model.dat$splits$valid@x,
    W_valid = model.dat$masks$valid,
    y = model.dat$depart,
    trace = 2,
    quiet = F,
    seed = 2023,
    early.stopping = 1,
    lambda.beta.grid = seq(10,0,length.out=20),
    max_cores = 5
  ) -> fit_bixi2
  
  fit_bixi2$hparams
  fit_bixi2$fit$X = X
  fit_bixi2$fit$beta[,1:5]
  print_performance(model.dat, fit_bixi2$fit, error_metric$rmse, F, "SoftImpute", F, 3)
  apply(fit_bixi2$fit$beta, 1, summary) |> as.data.frame() |> 
    setNames(colnames(X)) |>  round(2) |>
    t() |> as.data.frame() |>  
    mutate(Prop_non_zero = rowSums(round(fit_bixi2$fit$beta,6) > 0) / ncol(fit_bixi2$fit$beta)) |>
    mutate(all_zeros = Prop_non_zero == 0) |> 
    arrange(desc(Prop_non_zero)) |> 
    round(2) |> 
    arrange(desc(Median)) |>  kable()
  
  Prop_non_zero = rowSums(round(fit_bixi2$fit$beta,6) > 0) / ncol(fit_bixi2$fit$beta)
  sig_cols <- colnames(X)[Prop_non_zero>.2]
  
  fit_bixi2$fit$beta |> 
    as.data.frame() |> 
    setNames(as.character(model.dat$dates)) |> 
    t() |> 
    as.data.frame() |> 
    setNames(colnames(X)) |> 
    mutate(Time = model.dat$dates) |> 
    pivot_longer(-Time, names_to="Covariate", values_to="Value") |> 
    filter(Covariate %in% sig_cols) |> 
    ggplot(aes(Time, Value, color=Covariate)) +
    geom_point() +
    facet_wrap(~Covariate, scales = "free_y") +
    labs(
      title = "Station Covariate Coefficients",
      x = "Time",
      y = "Beta Coefficient"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "none"
    ) + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b")
    
  
    model.dat$Z |> 
      scale() |> 
      as.data.frame() |> 
      select(-holiday) |> 
      mutate(Time = model.dat$dates) |> 
      pivot_longer(-Time, names_to="Covariate", values_to="Value") |>  
      ggplot(aes(Time, Value, color=Covariate)) +
      geom_point() +
      facet_wrap(~Covariate, scales = "free_y") +
      labs(
        title = "Time Covariates",
        x = "Time",
        y = "Covariate Value"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(face = "bold"),
        legend.position = "none"
      ) + 
      scale_x_date(date_breaks = "1 month", date_labels = "%b")
    
    
    fit_bixi2$fit$M = unsvd(fit_bixi2$fit)
    M_melted <- melt(fit_bixi2$fit$M)
    M_melted$Date <- rep(model.dat$dates, each = nrow(model.dat$depart))
    
    # Create a heatmap using ggplot2
    ggplot(data = M_melted, aes(x = Var1, y = Date, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                           midpoint = 0, limit = c(min(M_melted$value), max(M_melted$value)), 
                           space = "Lab", name="Value") +
      theme_minimal() +
      labs(title = "Heatmap of Low-Rank Matrix M", x = "Station", y = "Date") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_flip() +
      scale_y_date(date_breaks = "7 days", date_labels = "%d/%m")
    
    
    data.frame(Y = as.vector(model.dat$splits$test[model.dat$masks$test==0]),
               M = as.vector(fit_bixi2$fit$M[model.dat$masks$test==0])) |> 
      ggplot(aes(Y, M)) +
      geom_point(alpha = 0.5) +
      # geom_smooth(method = "lm", col = "red") +
      geom_abline(intercept = 0, slope = 1, color = "red") +  
      labs(title = "Scatter Plot of Vectorized M against Vectorized Y of the test set",
           x = "M (vectorized)",
           y = "Y (vectorized)") +
      theme_minimal()
    
    data.frame(Y = as.vector(model.dat$splits$test[model.dat$masks$test==0]),
               M = as.vector(fit_bixi2$fit$M[model.dat$masks$test==0])) |>
      mutate(Xbeta = (X %*% fit_bixi2$fit$beta)[model.dat$masks$test==0] ) |> 
      mutate(preds = Xbeta + M) |> 
      ggplot(aes(Y, preds)) +
      geom_point(alpha = 0.5) +
      # geom_smooth(method = "lm", col = "red") +
      geom_abline(intercept = 0, slope = 1, color = "red") +  
      labs(title = "Scatter Plot of Vectorized estimates against Vectorized Y of the test set",
           x = "estimates (vectorized)",
           y = "Y (vectorized)") +
      theme_minimal()
    
    data.frame(Y = as.vector(model.dat$splits$test[model.dat$masks$test==0])) |>
      mutate(Xbeta = (X %*% fit_bixi2$fit$beta)[model.dat$masks$test==0] ) |> 
      ggplot(aes(1:length(Xbeta), Xbeta)) +
      geom_point(alpha = 0.5) +
      # geom_smooth(method = "lm", col = "red") +
      #geom_abline(intercept = 0, slope = 1, color = "red") +  
      labs(title = "X beta values for the test set arranged by station and date",
           x = "",
           y = "X * beta (vectorized)") +
      theme_minimal()
    
    
  #------------------------------------------------
#   results <- CASMC_var_selection(
#     y_train = model.dat$splits$train,
#     y_valid = model.dat$splits$valid,
#     Y = model.dat$depart,
#     X = X,
#     W_valid = model.dat$masks$valid,
#     track = F,
#     return_best_model = TRUE,
#   )
#   
#   results$res |>
#     select(-variables) |>
#     arrange(validation_error) |>
#     kable() |>
#     print()
#   fiti <- results$fit$fit
#   fiti$X <- results$fit$X
#   
#   all_res <- rbind(all_res,
#                    print_performance(
#                      model.dat,
#                      fiti,
#                      error_metric$rmse,
#                      F,
#                      paste0("CASMC(", case, ")"),
#                      F,
#                      3
#                    ))
#   print(all_res)
# #}



simpute_fit <- simpute.cv(
  Y_train = as.matrix(model.dat$splits$train),
  y_valid = model.dat$splits$valid@x,
  W_valid = model.dat$masks$valid,
  y = as.matrix(model.dat$depart),
  n.lambda = 20,
  trace = FALSE,
  print.best = TRUE,
  tol = 5,
  thresh = 1e-6,
  rank.init = 2,
  rank.limit = 30,
  rank.step = 2,
  maxit = 800,
  seed= 2023
  
)
print_performance(model.dat, simpute_fit, error_metric$rmse, TRUE, "SoftImpute", F, 3)

# 
# all_res  |>
#   rbind(print_performance(model.dat, simpute_fit, error_metric$rmse, TRUE, "SoftImpute", F, 3)) |> 
#   #arrange(XBeta) |>
#   #mutate(XBeta = paste0("[", 1:nrow(all_res), "]", XBeta)) |>
#   #arrange(M) |>
#   #mutate(M = paste0("[", 1:nrow(all_res), "]", M)) |>
#   #arrange(Y) |>
#   #mutate(Y = paste0("[", 1:nrow(all_res), "]", Y)) |>
#   dplyr::select(-Y, -M, -XBeta) |>
#   arrange(train) |>
#   mutate(train = paste0("[", 1:nrow(all_res), "]", train)) |>
#   arrange(test) |>
#   mutate(test = paste0("[", 1:nrow(all_res), "]", test))  ->
#   all_res
# 
# 
# kable(all_res, format = "simple") |>  print()





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #-----------------------------------------------------------------------------
# aresults <- list()
# i = b = 1
# 
# X_r <- reduced_hat_decomp(model.dat$X, 1)
# X_r$rank
# 
# sparm = list(NULL,NULL,"")
# for(sparm in list(list(NULL,NULL,""),
#                   list(model.dat$similarity.A, NULL, "_with_time_laplacian")
#                   #list(NULL, model.dat$sim_col),
#                   #list(model.dat$sim_row, model.dat$sim_col)
#                   )
#     ){
# 
#   X <- model.dat$X[,c(1), drop=FALSE]|> scale()
#   #X = cbind(X[,1:4]^2)
#   #------------------------
#   X <- cbind(
#     #X,
#     scale(X[,1:3])^2,
#     log(X[,1, drop=FALSE]+abs(min(X[,1]))+1)
#     #matrix(X[,4]>0,ncol= 1)
#   )  
#   
#   #---------------------
#   cor(X)
#   # for(b in 1:5){
#   model.dat <- load_bixi_dat(transpose = F, scale_response = F, scale_covariates = F,
#                              testp = 0.5, validp = 0.2, seed=b)$model
#   model.dat$X <- model.dat$X |> scalers("minmax")
#   model.dat$X <- remove_collinear_cols(model.dat$X, 0.7)
#   model.dat$X <- reduced_hat_decomp(model.dat$X, .01, 0.98)$X
#   #X <-  model.dat$X[,-c(1,5,9,12,16)] |> scalers("minmax")
#   #X <- X[,-c(5,7,9)]
#   #X_r = reduced_hat_decomp(X,.01, 0.99)
#   #X = X[,1:5]
#   #X[1,]
#   #dim(X_r$X)
#   #X_r$rank
#   #X <- model.dat$X[,c(1,2, 3,4,6,7,8)]
#   
#   results <- CASMC_var_selection(
#     y_train = model.dat$splits$train,
#     y_valid = model.dat$splits$valid,
#     Y = model.dat$depart,
#     X = model.dat$X,
#     W_valid = model.dat$masks$valid,
#     track = TRUE
#   )
#   
#   results |> 
#     select(-variables) |> 
#     arrange(validation_error)
#   
#   #-----------------------------------------------------------------------------
#   start_time = Sys.time()
#   best_fit = CASMC_cv_rank(
#     y_train = model.dat$splits$train,
#     X = X,#[,1:5, drop=FALSE],
#     
#     y_valid = model.dat$splits$valid@x,
#     W_valid = model.dat$masks$valid ,
#     #y = model.dat$depart,
#     trace = F,
#     max_cores = 30,
#     thresh = 1e-6,
#     lambda.a = 0.01,
#     S.a = sparm[[1]],
#     lambda.b = 0.2,
#     S.b = sparm[[2]],
#     #n.lambda = n.lambda,
#     #rank.limit = rank.limit,
#     maxit = 200,
#     r_min = 0,
#     rank.init = 2,
#     rank.step = 4,
#     print.best = TRUE,
#     seed = 2023,
#     track  = T
#   )
#   best_fit$rank_M
#   test_error <- error_metric$rmse
#   fit1 = best_fit$fit
#   sout = best_fit
#   # get estimates and validate
#   sout$M = unsvd(fit1)
#   sout$beta =  fit1$beta
#   apply(sout$beta, 1, summary) |> print()
#   
#   sout$estimates = sout$M + X %*% (sout$beta)
#   
#   hist(sout$beta[1,])
#   
#   plot(1:length(sout$estimates), as.vector(sout$estimates))
#   
#   #plot(1:ncol(sout$beta), sout$beta[3,])
#   
#   if(b > 1) old_results = results
#   results = list(model = paste0("CASMC_rank_all",sparm[[3]]))
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
  #------------------------------
#=========================================================================================
#   
# X <- X[,1:3]
#   for(b in 1:5){
#     model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
# start_time = Sys.time()
# 
# best_fit = CASMC_cv_rank(
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
# results = list(model = paste0("CASMC_rank_1:3_",sparm[[3]]))
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
# best_fit = CASMC_cv_rank(
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
# results = list(model = paste0("CASMC_rank_2_",sparm[[3]]))
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
# best_fit = CASMC_cv_L2(
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
# results = list(model = paste0("CASMC_rank_L2_all_",sparm[[3]]))
# results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
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
# 
