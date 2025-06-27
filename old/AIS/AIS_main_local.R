#setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/code_files/AIS")
setwd("E:/HEC/Youssif/github_folder_DONT_PUSH/HEC_MAO_COOP/code_files/AIS")



source("./AIS_fit_mine_local.R")
source("./AIS_fit_mine_cov_local.R")
source("./Approximate_SVT_local.R")
source("./power_method_local.R")
source("../Mao_import_lib_local.R")



gen.dat <- generate_simulation_data_ysf(2,600,600,10,10, missing_prob = 0.9,coll=TRUE)
# W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
# Y_train = (gen.dat$Y * W_valid)
# Y_valid = gen.dat$Y[W_valid==0]
#-----------------------------------
start_time <- Sys.time()
fit_ = AIS_cov_fit(gen.dat$Y, gen.dat$X, lambda=8.3799, trace.it = F, tol=1e-3,
               maxIter = 100, decay = 0.8, maxR = 10)
print(paste("Execution time is",
            round(as.numeric(difftime(Sys.time(),
                                      start_time,units = "secs"))), "seconds"))

preds <- fit_$U0 %*% (fit_$S * fit_$V)


test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])


#----------------------------------
start_time <- Sys.time()
fit_ = AIS_fit(gen.dat$Y, lambda=12.422, trace.it = F, tol=1e-3,
               maxIter = 300, decay = 0.5, maxR = 15)
print(paste("Execution time is",
            round(as.numeric(difftime(Sys.time(),
                                      start_time,units = "secs"))), "seconds"))
preds <- fit_$U0 %*% (fit_$S * fit_$V)
test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])



Y = gen.dat$Y; lambda=12.422; maxR = NA; maxIter=100; tol=1e-6;
thresh=1e-3;decay=0.9; trace.it=FALSE; exact = FALSE


W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]
# fiting Alogirthm 1
lambda1 = 23.333
lambda2 = 8.379
max.rank = 16
beta_partial = solve(t(gen.dat$X) %*% gen.dat$X + 
                       diag(lambda1, ncol(gen.dat$X))) %*% t(gen.dat$X)
start_time <- Sys.time()
sout <- simpute.als.cov(gen.dat$Y, gen.dat$X, beta_partial,
                        J = max.rank, thresh =  1e-3, lambda= lambda2,
                        trace.it = F,warm.start = NULL)
print(paste("Execution time is",
            round(as.numeric(difftime(Sys.time(),
                          start_time,units = "secs"))), "seconds"))
sout$A_hat = sout$u %*% (sout$d * t(sout$v))
print(paste("Test error =", round(test_error(sout$A_hat[gen.dat$W==0],
                                             gen.dat$A[gen.dat$W==0]),5)))
#-------------------------------------------------------
start_time <- Sys.time()
sout <- simpute.cov.cv(Y_train, gen.dat$X, W_valid, Y_valid, 
                       trace=TRUE, rank.limit = 30, lambda1=0,n1n2 = 1, 
                       warm=NULL,tol = 2, print.best = FALSE)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), 
                            start_time,units = "secs"))), "seconds"))
print(paste("Test error =", round(test_error(sout$A_hat[gen.dat$W==0],
                                             gen.dat$A[gen.dat$W==0]),5)))
print(sprintf("Optimal hyperparameters are: lambda2 = %.4f, r=%d", sout$lambda,
              sout$rank.max))