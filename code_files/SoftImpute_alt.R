
source("./code_files/Mao_import_lib.R")
library(softImpute)
dim = c(800)
missingness = 0
i=1
coll=TRUE

if(missingness == 0){
   gen.dat <- generate_simulation_data_mao(n1=dim[i],n2=dim[i],m=5,r=10, seed=2023)
}else
   gen.dat <- generate_simulation_data_ysf(2,dim[i],dim[i],10,10, missing_prob = missingness,coll=coll)

lambda.1_grid = seq(0,3,length=20)
lambda.2_grid = seq(.9, 0, length=20)
alpha_grid = seq(0.992, 1, length=10)

cv.out <- Mao.cv(gen.dat$A, gen.dat$X, gen.dat$Y, gen.dat$W,
                 n_folds=5, 
                 lambda.1_grid = lambda.1_grid,
                 lambda.2_grid = lambda.2_grid,
                 alpha_grid = alpha_grid,
                 numCores = 1,n1n2_optimized = TRUE,theta_estimator = theta_default)
mao.out <- Mao.fit(gen.dat$Y, gen.dat$X, gen.dat$W, cv.out$best_parameters$lambda.1, 
                   cv.out$best_parameters$lambda.2, cv.out$best_parameters$alpha, 
                   theta_estimator = theta_default)

test_error(mao.out$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
#-----------------------------------------------------------
# 1. Soft Impute with/without X
with_X = FALSE

if (with_X==TRUE){new_Y = gen.dat$Y - gen.dat$X %*% mao.out$beta_hat}else new_Y = gen.dat$Y
new_Y[gen.dat$W==0] = NA

#xs <- as(new_Y, "Incomplete")
lam0 <- lambda0(new_Y)
lamseq <- exp(seq(from=log(lam0), to=log(1), length=20))

fits <- as.list(lamseq)
ranks <- as.integer(lamseq)
rank.max <- 10
warm <- NULL


for(i in seq(along=lamseq)) {
   fiti <- softImpute(new_Y, lambda=lamseq[i], rank=rank.max, warm=warm)
   ranks[i] <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
   rank.max <- min(ranks[i]+2, 50)
   if (with_X==TRUE){ soft_estim  = complete(new_Y, fiti) +  gen.dat$X %*% mao.out$beta_hat}else soft_estim  =  complete(new_Y, fiti)
   
   err = test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   warm <- fiti # warm start for next 
   fits[[i]] <- fiti
   cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
               i, lamseq[i], rank.max, ranks[i], err))
}
#---------------------------------------------
fitslo <- simpute.svd.cov(new_Y, gen.dat$X, gen.dat$W,lambda = 7, trace.it=FALSE, warm.start = fitslo)

v=as.matrix(fitslo$v)
vd=v*outer(rep(1,nrow(v)),fitslo$d)
soft_estim = fitslo$u %*% t(vd)  + gen.dat$X %*% fitslo$beta.estim
test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])

fitslo$d
#-----------------------------------------------------------------------------------
# with covariates - loop
new_Y = gen.dat$Y
new_Y[gen.dat$W==0] = NA

#xs <- as(new_Y, "Incomplete")
lam0 <- lambda0(new_Y)
lam0 <- lambda0.cov(new_Y, gen.dat$X, gen.dat$W)
lam0 <- 40
lamseq <- seq(from=10, to=0, length=20)

fits <- as.list(lamseq)
ranks <- as.integer(lamseq)
rank.max <- 10
warm <- NULL
best_fit <- list(error=Inf, rank=NA, lambda=NA, rank.max=NA)

for(i in seq(along=lamseq)) {
   fiti <- simpute.svd.cov(new_Y, gen.dat$X, gen.dat$W,lambda = lamseq[i], rank=rank.max, warm.start = warm)
   ranks[i] <- sum(round(fiti$d, 4) > 0) # number of positive sing.values
   rank.max <- min(ranks[i]+2, 50)
   # get estimates
   v=as.matrix(fiti$v)
   vd=v*outer(rep(1,nrow(v)),fiti$d)
   soft_estim = fiti$u %*% t(vd)  + gen.dat$X %*% fiti$beta.estim
   #----------------------------
   err = test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])
   warm <- fiti # warm start for next 
   fits[[i]] <- fiti
   cat(sprintf("%2d lambda=%9.5g, rank.max = %d  ==> rank = %d, error = %.5f\n",
               i, lamseq[i], rank.max, ranks[i], err))
   #-------------------------
   # register best fir
   if(err < best_fit$error){
      best_fit$error = err
      best_fit$rank_B = ranks[i]
      best_fit$rank_A = qr(soft_estim)$rank
      best_fit$lambda = lamseq[i]
      best_fit$rank.max = rank.max
   } 
}
print(best_fit)
#-------------------------------



fits <- softImpute(new_Y, trace=FALSE)
full.Y <- complete(new_Y, fits)
soft_estim <- full.Y + gen.dat$X %*% mao.out$beta_hat
test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0])


dim(new_Y)
X <- gen.dat$X

qr(full.Y)$rank
fits$d
Y <- gen.dat$Y

