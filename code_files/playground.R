
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

gen.dat <- generate_simulation_data_ysf(1,500,500,5,10, missing_prob = 0.9,coll=T,seed=3023)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_valid <- gen.dat$Y[W_valid==0]

Y_train <- gen.dat$Y
Y_train[Y_train==0] = NA
y = Y_train %>% as.matrix()
ys <- as(y, "Incomplete")


fits <- simpute.als.fit_Incomplete(ys, gen.dat$X, trace=TRUE, lambda=3)
fits$rank
fits$lambda
system.time({fitss <- simpute.als.fit_Incomplete(ys, J = 3, lambda=1.9)})







out <- coop_find_rho(gen.dat, W_valid,  print_best = TRUE,early_maxiter = 50,max_cores = 10,
                     rho.grid = seq(0.1,0.99,length.out=10))
out$time_in_minutes
out$rho

test_error(out$fit$preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])


out <- coop_fit(gen.dat$Y, gen.dat$X, gen.dat$Z, gen.dat$W, W_valid,
                Y_valid, rho=0.1, tol=2,trace_fin = F, verbose=FALSE, early_stopping = TRUE,
                patience=5, maxiter = 15, seed = 3023,print.best = T)

test_error(out$preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
test_error(out$preds[W_valid==0], gen.dat$A[W_valid==0])

out$valid_error

require(softImpute)
require(Matrix)
Y_train <- gen.dat$Y
Y_train[Y_train==0] = NA
fits <- softImpute(Y_train, rank.max=3, lambda=1.9, trace=TRUE)
y = Y_train %>% as.matrix()
ys <- as(y, "Incomplete")
ysc <- biScale(ys, col.scale=FALSE,row.scale=FALSE,trace=TRUE)
system.time({fitss <- softImpute(ysc, rank.max = 3, lambda=1.9)})
system.time({fits <- softImpute(Y_train, rank.max = 3, lambda=1.9)})

system.time({complete(ysc, fitss)})
fitss <- deBias(ysc, fitss)
preds <- fitss$u %*% (fitss$d * t(fitss$v))
test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])

system.time({complete(Y_train, fits)})
fits <- deBias(Y_train, fits)
preds <- fits$u %*% (fits$d * t(fits$v))
test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
qr(preds)$rank

#=====================
i=row(y)[!is.na(y)]
j=col(y)[!is.na(y)]
value=y[!is.na(y)]
cbind(i,j,value)



## -----------------------------------------------------------------------------
Incomplete(i=i,j=j,x=value)

