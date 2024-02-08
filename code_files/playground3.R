
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")


gen.dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.9,coll=T)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]
X <- gen.dat$X
lambda2 = 31
max.rank = 15
beta_partial = solve(t(gen.dat$X) %*% gen.dat$X) %*% t(gen.dat$X)


# folds <- k_fold_cells(nrow(gen.dat$Y), ncol(gen.dat$Y), 3, gen.dat$W)[[1]]

X_reduced = reduced_hat_decomp(gen.dat$X, 1e-2)

# best_fit = simpute.cov.cv_splr_no_patience(Y_train, X_reduced$svdH, Y_valid, W_valid,warm = best_fit$best_fit,
#                                            trace = T, rank.limit=30,rank.step=4,patience = 3)


start_time = Sys.time()
best_fit = simpute.cov.Kf_splr_no_patience_v2(gen.dat$Y, X_reduced$svdH, gen.dat$W, n_folds=10,
                                            trace = F, rank.limit=30,rank.step=4,patience = 1)


yfill = gen.dat$Y
fits = best_fit$best_fit
best_fit$B_hat = fits$u %*% (fits$d * t(fits$v))
yfill[gen.dat$Y==0] = (best_fit$B_hat)[gen.dat$Y==0]
fits.out = list(u=fits$u, d=fits$d, v=fits$v, beta.estim=beta_partial %*% yfill)
X_svd = X_reduced$X
beta_partial = MASS::ginv(t(X_svd) %*% X_svd) %*% t(X_svd)


set.seed(2023);fits2 <- simpute.als.cov(gen.dat$Y, X_svd, beta_partial,J = best_fit$rank_B, thresh =  1e-6,
                                        lambda= best_fit$lambda,trace.it = F,warm.start = fits.out, maxit=100)

fits2$M = fits2$u %*% (fits2$d * t(fits2$v))
fits2$A_hat = fits2$M  + X %*% fits2$beta.estim
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))


test_error(fits2$beta.estim, gen.dat$beta)
print(paste("Test error =", round(test_error(fits2$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (fits2$A_hat[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
print(paste("Test error =", round(test_error(fits2$M, gen.dat$B),5)))
test_error(best_fit$B_hat, gen.dat$B)
#######################################################
# original
start_time = Sys.time()
sout <- simpute.cov.cv(Y_train, gen.dat$X, W_valid, Y_valid, trace=FALSE, rank.limit = 30, 
                       print.best=FALSE, rank.step=4, type="als", lambda1=0, tol=2)
sout <- simpute.cov.cv.lambda1(Y_train, gen.dat$X, W_valid, Y_valid, sout$lambda, sout$rank.max, print.best = FALSE,
                               trace=FALSE, lambda1.grid = seq(0,20,length.out=20) ,n1n2 = 1, warm=NULL)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(sout$beta_hat, gen.dat$beta)
test_error(sout$B_hat, gen.dat$B)
#################################################################

start_time = Sys.time()
soutk <- simpute.cov.kfold(gen.dat$Y, gen.dat$X, gen.dat$W, n_folds = 3, print.best = FALSE,
                         trace=FALSE, rank.limit = 30, lambda1=0,n1n2 = 1, warm=NULL,tol = 2)
soutk <- simpute.cov.kfold.lambda1(gen.dat$Y, gen.dat$X, gen.dat$W, soutk$lambda2, n_folds = 3, print.best = FALSE, 
                                  trace=FALSE,lambda1.grid = seq(0,20,length.out=20) ,n1n2 = 1, warm=NULL,
                                  J=c(soutk$J))
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
test_error(soutk$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(soutk$beta_hat, gen.dat$beta)
test_error(soutk$B_hat, gen.dat$B)
