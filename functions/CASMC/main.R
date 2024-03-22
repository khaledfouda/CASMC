# this file illustrates an example of using the methods in the folder

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data

gen.dat <-  generate_simulation_data_mao(800,800,10,10, "MAR", 2024)
gen.dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.9,coll=T)
X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
y = yfill = gen.dat$Y#Y_train
y[y==0] = NA
y <- as(y, "Incomplete")
xbeta.sparse = y
lambda2 =  15.66312
max.rank = 3

#-----------------------------------------------------------------------
# 1. Fit the model without cross validation:
# files: fit_covariates.R; fit_n_covariates_fixed_rank_beta.R

start_time <- Sys.time()
set.seed(2020);fits <- CASMC_fit(y=y, svdH=X_r$svdH,  trace=F, J=max.rank,
                                            thresh=1e-6, lambda=lambda2, init = "naive",
                                            final.svd = T,maxit = 500, warm.start = NULL)
xbeta.sparse@x = fits$xbeta.obs
# warm.start above expects u,d,v,xbeta.obs
fit4 = SZIRCI(xbeta.sparse, X_r$X,  X_r$rank, maxit=300,
                                trace.it = F,final.trim = F,thresh = 1e-6)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

# get estimates and validate
beta =  fit4$u %*% (fit4$d * t(fit4$v))
M = fits$u %*% (fits$d * t(fits$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] )
test_error(fits$xbeta.obs, (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta),5)))
test_error(M, gen.dat$B)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
#----------------------------------------------------
# second version option 1
start_time <- Sys.time()
set.seed(2020);fits2 <- CASMC_fit_v3(y=y, X=X_r$X, svdH=X_r$svdH,  trace=F, J=max.rank,
                                 thresh=1e-6, lambda=lambda2, init = "naive",
                                 final.svd = T,maxit = 500, warm.start = NULL)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

# get estimates and validate
beta =  fits2$beta

M = fits2$u %*% (fits2$d * t(fits2$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] )
test_error(fits2$xbeta.obs, (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta),5)))
test_error(M, gen.dat$B)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))



#------------------------------------------------------------------------------------
# [optional] prepare warm.start for the second fit function:

svdX = fast.svd(X_r$X)
Ux = svdX$u
Vx = svdX$d * t(svdX$v)
X0 = ginv(t(Vx)%*%Vx) %*% t(Vx)
warm.start = list()
warm.start$X1 = X0 %*% t(Ux)
warm.start$X2 = X0 %*% Vx
B = t( ginv(X_r$X) %*% naive_MC(as.matrix(xbeta.sparse))) # B = (X^-1 Y)'
warm.start$Bsvd = fast.svd(B)
fit4 = SZIRCI(xbeta.sparse, X_r$X,  X_r$rank, maxit=300, warm.start = warm.start,
                                 trace.it = T,final.trim = F,thresh = 1e-6)
#----------------------------------------------------------------------------------------------------------
# 2. cross-validation with train/test split
# file: split_cross_validation.R

# prepare the data
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]
# fit
start_time <- Sys.time()
best_fit = CASMC_cv_holdout(Y_train, X_r, Y_valid, W_valid, gen.dat$Y,
                               trace=F, thresh=1e-6,n.lambda = 30, rank.limit = 20)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fit1 = best_fit$fit1
fit2 = best_fit$fit2
# get estimates and validate
beta =  fit2$u %*% (fit2$d * t(fit2$v))
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta.x)[gen.dat$Y!=0] )
# test_error(fit1$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta.x),5)))
test_error(M, gen.dat$B)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
best_fit$rank.max
best_fit$lambda
#------------------------------------------------------------------------------------
# K-fold cross-validation
test_error <- RMSE_error
start_time <- Sys.time()
set.seed(2023)
best_fit2 = CASMC_cv_kfold(gen.dat$Y, X_r, gen.dat$W, trace=T,print.best = TRUE,
                                n.lambda = 20, n_folds = 3, rank.limit=30, rank.step=2)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fit1 = best_fit2$fit1
fit2 = best_fit2$fit2
# get estimates and validate
beta =  fit2$u %*% (fit2$d * t(fit2$v))
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + X_r$X %*% t(beta)

#test_error <- adjusted_unexplained_variance

test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta.x)[gen.dat$Y!=0] )
# test_error(fit1$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta.x),5)))
test_error(M, gen.dat$B)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
best_fit2$rank.max
best_fit2$lambda
#---------------------------------------------------
# the following also does cross validation but using the functions defined
# in compare_models

valid.dat = list()
valid.dat$W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
valid.dat$Y_train <- gen.dat$Y * valid.dat$W_valid
valid.dat$Y_valid <- gen.dat$Y[valid.dat$W_valid==0]
splr.dat = reduced_hat_decomp(gen.dat$X, 1e-2)
compare_softImpute_splr_Kfold(gen.dat, splr.dat, 3,2, 20,30) -> results
results = rbind(results, compare_softImpute_splr(gen.dat, valid.dat, splr.dat,
                                                 2, 20,30))
results
#-=-------------------------------------------------------------------------------------------
# this is a test of the second function 
# 
B <- 50
with_noise = TRUE
error1 <- error2 <- numeric(B)
time1 <- time2 <- 0
for(b in 1:B){
   
   if(b%%10==0) print(b)
   gen.dat <- generate_simulation_data_ysf(2,400,300,10,10, missing_prob = 0.9,coll=T,seed=b)
   X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
   Y = A = X_r$X %*% gen.dat$beta
   if(with_noise) Y <- Y + matrix(rnorm(dim(Y)[1]*dim(Y)[2]),dim(Y)[1],dim(Y)[2])
   Y[gen.dat$W == 0] = NA 
   Y <- as(Y, "Incomplete")
   
   start_time <- Sys.time()
   fit = simpute.als.splr.fit.beta(Y, X_r$X,  X_r$rank, maxit=300,
                                    trace.it = F,final.trim = F,thresh = 1e-6)
   
   beta.estim <- fit$u %*% (fit$d * t(fit$v))
   A.estim <- X_r$X %*% t(beta.estim)
   error1[b] <- RMSE_error(A.estim[gen.dat$W==0], A[gen.dat$W==0])
   
   time1 <- time1 + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
   
   start_time <- Sys.time()
   fit2 <- softImpute(Y, rank.max=X_r$rank, type="als", thresh=1e-6, trace.it=F,maxit = 300) 
   A.estim2 <- fit2$u %*% (fit2$d * t(fit2$v))
   error2[b] <- RMSE_error(A.estim2[gen.dat$W==0], A[gen.dat$W==0])
   
   time2 <- time2 + as.numeric(difftime(Sys.time(), start_time,units = "secs"))

}

print(paste0("The new model has on average (with sd) rmse of ",round(mean(error1),4),
             "(",round(sd(error1),4),") and took ", round(time1),
             " seconds to run ", B, " simulations"))
print(paste0("The SoftImpute model has on average (with sd) rmse of ",round(mean(error2),4),
             "(",round(sd(error2),4),") and took ", round(time2),
             " seconds to run ", B, " simulations"))


