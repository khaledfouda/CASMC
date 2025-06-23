setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
library(BKTR)
source("./code_files/import_lib.R")
source("./BIXI/data-raw/bixi_data.R")

#----
model.dat <- load_bixi_dat(transpose = T, scale_response = F)$model
#-----------------

#----------------------------------------------------------------------------------
# apply models:
model.dat$X[1,]
# CAMC
aresults <- list()
i = 1

X_r <- reduced_hat_decomp(model.dat$X, 1)
X_r$rank

sparm = list(NULL,NULL,"")
for(sparm in list(list(NULL,NULL,""),
                  list(model.dat$similarity.A, NULL, "_with_time_laplacian")
                  #list(NULL, model.dat$sim_col),
                  #list(model.dat$sim_row, model.dat$sim_col)
                  )
    ){

  X <- model.dat$X[,c(1), drop=FALSE]|> scale()
  #X = cbind(X[,1:4]^2)
  #------------------------
  X <- cbind(
    #X,
    scale(X[,1:3])^2,
    log(X[,1, drop=FALSE]+abs(min(X[,1]))+1)
    #matrix(X[,4]>0,ncol= 1)
  )  
  
  #---------------------
  cor(X)
  for(b in 1:5){
  model.dat <- load_bixi_dat(transpose = F, scale_response = T, scale_covariates = F,
                             testp = 0.1, validp = 0.2, seed=b)$model
  
  X <-  model.dat$X[,-c(1,5,9,12,16)] |> scalers("minmax")
  X <- X[,-c(5,7,9)]
  X_r = reduced_hat_decomp(X,.01, 0.99)
  X = X[,1:5]
  X[1,]
  dim(X_r$X)
  X_r$rank
  #X <- model.dat$X[,c(1,2, 3,4,6,7,8)]
  
  start_time = Sys.time()
  best_fit = CAMC_cv_rank(
    y_train = model.dat$splits$train,
    X = X,#[,1:5, drop=FALSE],
    
    y_valid = model.dat$splits$valid@x,
    W_valid = model.dat$masks$valid ,
    y = model.dat$depart,
    trace = F,
    max_cores = 30,
    thresh = 1e-6,
    lambda.a = 0.01,
    S.a = sparm[[1]],
    lambda.b = 0.2,
    S.b = sparm[[2]],
    #n.lambda = n.lambda,
    #rank.limit = rank.limit,
    maxit = 200,
    r_min = 0,
    rank.init = 2,
    rank.step = 4,
    print.best = TRUE,
    seed = 2023,
    track  = T
  )
  best_fit$rank_M
  test_error <- error_metric$rmse
  fit1 = best_fit$fit
  sout = best_fit
  # get estimates and validate
  sout$M = unsvd(fit1)
  sout$beta =  fit1$beta
  apply(sout$beta, 1, summary) |> print()
  
  sout$estimates = sout$M + X %*% (sout$beta)
  
  hist(sout$beta[1,])
  
  plot(1:length(sout$estimates), as.vector(sout$estimates))
  
  #plot(1:ncol(sout$beta), sout$beta[3,])
  
  if(b > 1) old_results = results
  results = list(model = paste0("CAMC_rank_all",sparm[[3]]))
  results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
  #results$lambda.1 = NA#sout$lambda.beta |> round(3)
  #results$lambda.2 = sout$lambda |> round(3)
  results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
                                  model.dat$splits$test@x) |> round(5)
  results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                   model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
    round(5)
  results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
                                   model.dat$splits$valid@x) |> round(5)
  
  results$rank = qr(sout$estimates)$rank
  if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
  }
  results[-1] = mapply(function(x)x/5, results[-1])
  
  
  aresults[[i]] <- results
  i = i +1
  #------------------------------
  
  
X <- X[,1:3]
  for(b in 1:5){
    model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
start_time = Sys.time()

best_fit = CAMC_cv_rank(
 y_train = model.dat$splits$train,
 X = X,
 y_valid = model.dat$splits$valid@x,
 W_valid = model.dat$masks$valid ,
 y = model.dat$depart,
 trace = F,
 max_cores = 30,
 thresh = 1e-6,
 lambda.a = 0.01,
 S.a = sparm[[1]],
 lambda.b = 0.2,
 S.b = sparm[[2]],
 #n.lambda = n.lambda,
 #rank.limit = rank.limit,
 maxit = 200,
 #rank.step = rank.step,
 print.best = TRUE,
 seed = 2023,
 track  = T
)
test_error <- error_metric$rmse
fit1 = best_fit$fit
sout = best_fit
# get estimates and validate
sout$M = unsvd(fit1)
sout$beta =  fit1$beta
apply(sout$beta, 1, summary) |> print()
sout$estimates = sout$M + X %*% (sout$beta)
if(b > 1) old_results = results

results = list(model = paste0("CAMC_rank_1:3_",sparm[[3]]))
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#results$lambda.1 = NA#sout$lambda.beta |> round(3)
#results$lambda.2 = sout$lambda |> round(3)
results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
                                model.dat$splits$test@x) |> round(5)
results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
 round(5)
results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
                                model.dat$splits$valid@x) |> round(5)

results$rank = qr(sout$estimates)$rank
if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
  }
  results[-1] = mapply(function(x)x/5, results[-1])
  

aresults[[i]] <- results
i = i +1
#-----------------------------

X <- X[,2, drop=FALSE]

for(b in 1:5){
  model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
start_time = Sys.time()

best_fit = CAMC_cv_rank(
  y_train = model.dat$splits$train,
  X = X,
  y_valid = model.dat$splits$valid@x,
  W_valid = model.dat$masks$valid ,
  y = model.dat$depart,
  trace = F,
  max_cores = 30,
  thresh = 1e-6,
  lambda.a = 0.01,
  S.a = sparm[[1]],
  lambda.b = 0.2,
  S.b = sparm[[2]],
  #n.lambda = n.lambda,
  #rank.limit = rank.limit,
  maxit = 200,
  #rank.step = rank.step,
  print.best = TRUE,
  seed = 2023,
  track  = T
)
test_error <- error_metric$rmse
fit1 = best_fit$fit
sout = best_fit
# get estimates and validate
sout$M = unsvd(fit1)
sout$beta =  fit1$beta
sout$estimates = sout$M + X %*% (sout$beta)
if(b > 1) old_results = results

results = list(model = paste0("CAMC_rank_2_",sparm[[3]]))
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#results$lambda.1 = NA#sout$lambda.beta |> round(3)
#results$lambda.2 = sout$lambda |> round(3)
results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
                                model.dat$splits$test@x) |> round(5)
results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                 model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
  round(5)
results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
                                 model.dat$splits$valid@x) |> round(5)

results$rank = qr(sout$estimates)$rank
if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
}
results[-1] = mapply(function(x)x/5, results[-1])


aresults[[i]] <- results
i = i +1
#------------------------------
X <- model.dat$X |> scale()
X = cbind(X[,1:4]^2)

for(b in 1:5){
  model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
start_time = Sys.time()

best_fit = CAMC_cv_L2(
  y_train = model.dat$splits$train,
  X = X,
  y_valid = model.dat$splits$valid@x,
  W_valid = model.dat$masks$valid ,
  y = model.dat$depart,
  trace = F,
  max_cores = 30,
  thresh = 1e-6,
  lambda.a = 0.01,
  S.a = sparm[[1]],
  lambda.b = 0.2,
  S.b = sparm[[2]],
  #n.lambda = n.lambda,
  #rank.limit = rank.limit,
  maxit = 200,
  #rank.step = rank.step,
  print.best = TRUE,
  seed = 2023,
  track = T
)
test_error <- error_metric$rmse
fit1 = best_fit$fit
sout = best_fit
# get estimates and validate
sout$M = unsvd(fit1)
sout$beta =  fit1$beta
apply(sout$beta, 1, summary) |> print()
sout$estimates = sout$M + X %*% (sout$beta)
if(b > 1) old_results = results

results = list(model = paste0("CAMC_rank_L2_all_",sparm[[3]]))
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#results$lambda.1 = NA#sout$lambda.beta |> round(3)
#results$lambda.2 = sout$lambda |> round(3)
results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
                                model.dat$splits$test@x) |> round(5)
results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                 model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
  round(5)
results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
                                 model.dat$splits$valid@x) |> round(5)

results$rank = qr(sout$estimates)$rank
if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
}
results[-1] = mapply(function(x)x/5, results[-1])


aresults[[i]] <- results
i = i +1

}
#-------------------------------------------------------------------------------
# Soft Impute
for(b in 1:5){
  model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
start_time = Sys.time()
sout <- simpute.cv(
 Y_train = as.matrix(model.dat$splits$train),
 y_valid = model.dat$splits$valid@x,
 W_valid = model.dat$masks$valid,
 y = as.matrix(model.dat$depart),
 trace = T,
 rank.limit = 30,
 tol = 10,
 print.best = FALSE,
 rank.step = 2,
 maxit = 700
)

if(b > 1) old_results = results

results = list(model = "SoftImpute")
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#results$lambda.1 = NA
#results$lambda.2 = sout$lambda |> round(3)
results$error.test = test_error(sout$estimates[model.dat$masks$test == 0],
                                model.dat$splits$test@x) |> round(5)
results$error.train = test_error(sout$estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                 model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
 round(5)
results$error.valid = test_error(sout$estimates[model.dat$masks$valid == 0],
                                 model.dat$splits$valid@x) |> round(5)

results$rank = sout$rank_M
if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
}
results[-1] = mapply(function(x)x/5, results[-1])

aresults[[i]] <- results
i = i +1
#-----------------------------------------------------------------------------------
for(b in 1:5){
  model.dat <- load_bixi_dat(transpose = T, scale_response = F, seed=b)$model
start_time = Sys.time()
estimates = naive_MC(as.matrix(model.dat$splits$train))
if(b > 1) old_results = results

results = list(model = "Naive")
results$time = round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
#results$lambda.1 = NA
#results$lambda.2 = NA
results$error.test = test_error(estimates[model.dat$masks$test == 0],
                                model.dat$splits$test@x) |> round(5)
results$error.train = test_error(estimates[model.dat$masks$test == 1 & model.dat$masks$obs == 1],
                                 model.dat$depart[model.dat$masks$test == 1& model.dat$masks$obs == 1]) |> 
 round(5)
results$error.valid = test_error(estimates[model.dat$masks$valid == 0],
                                 model.dat$splits$valid@x) |> round(5)
results$rank = qr(estimates)$rank
if(b > 1)
  results[-1] = mapply(sum, old_results[-1], results[-1])
}
results[-1] = mapply(function(x)x/5, results[-1])

aresults[[i]] <- results
i = i +1
#---------------------------------------------------------------------------------
do.call(rbind, lapply(aresults, function(x) data.frame(t(unlist(x))))) |> kable()


