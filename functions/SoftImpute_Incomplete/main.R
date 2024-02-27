# this file illustrates an example of using the methods in the folder

setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

# 0. prepare the data

gen.dat <- generate_simulation_data_ysf(2,800,700,10,10, missing_prob = 0.5,coll=F)
X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
y = yfill = gen.dat$Y#Y_train
y[y==0] = NA
y <- as(y, "Incomplete")
xbeta.sparse = y
lambda2 = 11.08761
max.rank = 5

#-----------------------------------------------------------------------
# 1. Fit the model without cross validation:
# files: fit_covariates.R; fit_n_covariates_fixed_rank_beta.R

start_time <- Sys.time()
set.seed(2020);fits <- simpute.als.fit_splr(y=y, svdH=X_r$svdH,  trace=F, J=max.rank,
                                            thresh=1e-6, lambda=lambda2, init = "naive",
                                            final.svd = T,maxit = 500, warm.start = NULL)
xbeta.sparse@x = fits$xbeta.obs
#warm.start above expects u,d,v,xbeta.obs
fit4 = simpute.als.splr.fit.beta(xbeta.sparse, X_r$X,  X_r$rank, maxit=300,
                                 trace.it = F,final.trim = F,thresh = 1e-6)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

# get estimates and validate
beta =  fit4$u %*% (fit4$d * t(fit4$v))
M = fits$u %*% (fits$d * t(fits$v))
A = M + X_r$X %*% t(beta)
test_error((X_r$X %*% t(beta))[gen.dat$Y!=0], (gen.dat$X %*% gen.dat$beta.x)[gen.dat$Y!=0] )
test_error(fits$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[gen.dat$Y!=0] )
print(paste("Test error =", round(test_error(t(beta), gen.dat$beta.x),5)))
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
fit4 = simpute.als.splr.fit.beta(xbeta.sparse, X_r$X,  X_r$rank, maxit=300, warm.start = warm.start,
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
best_fit = simpute.cov.cv_splr(Y_train, X_r, Y_valid, W_valid, y, trace=FALSE, thresh=1e-6)
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

start_time <- Sys.time()
best_fit2 = simpute.cov.Kf_splr(gen.dat$Y, X_r, gen.dat$W, trace=TRUE,print.best = TRUE)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fit1 = best_fit2$fit1
fit2 = best_fit2$fit2
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