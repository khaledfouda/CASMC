# this file illustrates an example of using the methods in the folder

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data

gen.dat <-  generate_simulation_data_mao(400,400,10,10, "MAR", 2024, cov_eff = F)

gen.dat <- generate_simulation_data_ysf(2,900,800,10,10, missing_prob = 0.8,coll=T,cov_eff = T,
                                        informative_cov_prop = 0.8,
                                        seed = 2024)


X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
y = yfill = gen.dat$Y#Y_train
y[y==0] = NA
y <- as(y, "Incomplete")
xbeta.sparse = y
lambda2 =  15.66312
max.rank = 3

X_r$X <- gen.dat$X
#-----------------------------------------------------------------------
# 1. Fit the model without cross validation:
similarity.a = generate_similarity_matrix(900, 2023)
similarity.b = generate_similarity_matrix(800, 2023)

ll <- computeLaplacian(similarity.a)
ll[ll==0] <- NA
ll2 <- as(ll, "Incomplete")
length(ll2@x)
length(ll)
suvC(gen.dat$X, t(gen.dat$beta), y@i, y@p)[1:5]

(gen.dat$X  %*% gen.dat$beta)[gen.dat$Y!=0][1:5]
bb = svd(t(gen.dat$beta))
suvC(gen.dat$X %*% bb$v, t(bb$d * t(bb$u)), y@i, y@p)[1:5]


start_time <- Sys.time()
set.seed(2020);fits2 <- CASMC_fit(y=y, X=X_r$X, svdH=X_r$svdH,  trace.it=F, J=max.rank, r = 10,
                                   thresh=1e-6, lambda=lambda2, lambda.beta=0.1,
                                   final.svd = T,maxit = 500, warm.start = NULL,
                                   lambda.a = 0.0, S.a=similarity.a,lambda.b = 2.0, S.b=similarity.b
                                  )
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

# fits2$Beta$d
# fits2$Beta$v[1:3,]
# 
# 
# fits2$Beta


# get estimates and validate
beta =  fits2$Beta$u %*% (fits2$Beta$d * t(fits2$Beta$v))

M = fits2$u %*% (fits2$d * t(fits2$v))
A = M + X_r$X %*% t(beta)
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta),5)))
test_error(M, gen.dat$M)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$O[gen.dat$W==0]),5)))
fits2$n_iter

#----------------------------------------------------------------------------------------------------------
# 2. cross-validation with train/test split
# file: split_cross_validation.R

# prepare the data
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_train[Y_train==0] = NA
Y_train <- as(Y_train, "Incomplete")
Y_valid = gen.dat$Y[W_valid==0]
# fit
start_time <- Sys.time()
best_fit = CASMC_cv_holdout(Y_train, X_r, Y_valid, W_valid, r = NULL, y=y,
                               trace=T, thresh=1e-6,n.lambda = 30, rank.limit = 20)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fit1 = best_fit$fit
# get estimates and validate
beta = as.matrix(fit1$Beta$u %*% (fit1$Beta$d * t(fit1$Beta$v)) )
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] )
# test_error(fit1$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta),5)))
test_error(M, gen.dat$M)
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$O[gen.dat$W==0]),5)))
best_fit$rank.max
best_fit$lambda
#---------------------------------------------------------------------------------
# ---- <r>
# prepare the data


# fit
start_time <- Sys.time()
best_fit2 = CASMC_cv_holdout_with_reg(Y_train, X_r, Y_valid, W_valid,  y=y, 
                                      lambda.beta.grid =  seq(0,3,length.out=30), max_cores = 30,
                                    trace=F, thresh=1e-6,n.lambda = 30, rank.limit = 20, track_beta = T) 
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))


#or
start_time <- Sys.time()
best_fit2 = CASMC_cv_holdout_with_r(Y_train, X_r, Y_valid, W_valid, r_min=0,  y=y,
                            trace=F, thresh=1e-6,n.lambda = 30, rank.limit = 20, track_r = T) 
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

print(best_fit2$r)
fit1 = best_fit2$fit
# get estimates and validate
beta =  as.matrix(fit1$Beta$u %*% (fit1$Beta$d * t(fit1$Beta$v)) )
M = fit1$u %*% (fit1$d * t(fit1$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta)[gen.dat$Y!=0] ) 
# test_error(fit1$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] )
print(paste("Test error =", round(test_error(A[gen.dat$W==0], gen.dat$O[gen.dat$W==0]),5))) 
test_error(M, gen.dat$M)
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta),5)))
best_fit2$r
best_fit2$rank.max
best_fit2$lambda
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
#-----------------------------------------------------------------------------------
# Repeat with the new method
set.seed(2023)
start_time <- Sys.time()
best_fit2 = CASMC_cv_kfold_v2(gen.dat$Y, X_r, gen.dat$W, trace=T,print.best = TRUE,#thresh = 1e-7,maxit = 300,
                           n.lambda = 20, n_folds = 3, rank.limit=30, rank.step=2)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fit1 = best_fit2$fit
# get estimates and validate
beta =  fit1$beta
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
#----------------------------------------
# K-Fold with CAMAC
# set.seed(2023)
start_time <- Sys.time()
fitiC <- CAMC_cv_kfold(gen.dat$Y,
                      X_r$X,
                      gen.dat$W,
                      n_folds = 5,
                      trace = TRUE,
                      rank.limit = 30,
                      print.best = TRUE,
                      #lambda.1_grid = rep(0,10),#lambda.1_grid,
                      rank.step = 2,
                      n.lambda = 20)
sout <- fitiC$fit2
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))


test_error((X_r$X %*% (sout$beta_hat))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta.x)[gen.dat$Y!=0] )
# test_error(fit1$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] )
print(paste("Beta error =", round(test_error((sout$beta_hat), gen.dat$beta.x),5)))
test_error(sout$B_hat, gen.dat$B)
print(paste("Test error =", round(test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sout$lambda1
sout$lambda2
sout$rank_A
sout$J
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


