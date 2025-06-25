library(tidyverse)
library(dplyr)
library(magrittr)
library(microbenchmark)


m = 10
n = 15
J = 5
lambda.M = 1.2

Dsq <- rnorm(J)
y <- rbinom(m*n, 1, .2) %>% matrix(m,n) 
y[y==0] <- NA
y %<>% as("Incomplete")
U <- rnorm(m*J) %>% matrix(m,J)
V <- rnorm(n*J) %>% matrix(n,J)
Dstar = (Dsq / (Dsq + lambda.M))
VDsq = t(Dsq * t(V))


source("./code_files/import_lib.R")

dat <- generate_simulated_data(600, 700, 10, 12, .8, F,
                               prepare_for_fitting = T)
     
start_time <- Sys.time()
set.seed(2020); fits <- simpute.cv(
  as.matrix(dat$fit_data$train),
  dat$fit_data$valid,
  dat$fit_data$W_valid,
  dat$Y,
  trace = T
)
# rank 9; rank_m = 7; lambda=0;

utils$error_metric$rmse(fits$estimates[dat$fit_data$W_valid==0], 
                        dat$fit_data$valid ) # better than new.xbeta



Xs <- qr(dat$X)
R <- qr.R(Xs)
Xs <- qr.Q(Xs)

fit3 <- CAMC3_fit(
  dat$fit_data$train,
  Xs,#dat$X,
  J = 9,
  lambda.M = .0000001,
  lambda.beta = .1,
  trace.it=T
)


XbetaValid = Xs %*% (fit3$beta)
MValid = fit3$u %*% (fit3$d * t(fit3$v))
pred = (XbetaValid+MValid)[dat$fit_data$W_valid==0]
utils$error_metric$rmse(pred, dat$fit_data$valid)
#--------------------------------------------
# cv


hpar <- CASMC_Lasso_hparams
hpar$M$rank.init = 2
#hpar$M$lambda.init = .000001

for (lamb in seq(0, 29, length=20)){
  
fitcv <- CAMC3_cv_M(
  dat$fit_data$train,
  Xs,
  dat$fit_data$valid,
  dat$fit_data$W_valid,
  y = dat$fit_data$Y,
  lambda_beta = 0,
  trace = 0,
  hpar = hpar
); 

print(paste(lamb, ": ",sum(fitcv$fit$beta ==0) / length(fitcv$fit$beta)))
}

E <- dat$fit_data$train - fitcv$fit$u %*% (fitcv$fit$d * t(fitcv$fit$v)) -
  Xs %*% fitcv$fit$beta
XtE <- crossprod(Xs, E)
lambda_max =  
  max(XtE + fitcv$fit$beta)
lambda_max

#-----------------------------------

old_lambda = 29
for (lamb in seq(29, 1, length=20)){
  fit3 <- CAMC3_fit(
    dat$fit_data$train,
    Xs,#dat$X,
    J = 10,
    lambda.M = 1,
    lambda.beta = lamb,
    trace.it=F
  )
  beta_ratio = sum(fit3$beta ==0) / length(fit3$beta)
  #print(paste(lamb, ": ",beta_ratio))
  
  if(beta_ratio < 1){
    print("Reached smalled beta")
    lambda_beta_max = old_lambda
    break
  }
  old_lambda = lamb 
}

#---------------------------------


fitcv2 <- CAMC_Lasso_cv(
  dat$fit_data$train,
  Xs,
  dat$fit_data$valid,
  dat$fit_data$W_valid,
  y = dat$fit_data$Y,
  #lambda_beta = 0.1,
  trace = 5,
  max_cores = 1,
  hpar = hpar
)


getdef
#-------------------------------------------

X <- Xs
y <- dat$fit_data$train
irow <- y@i
pcol <- y@p
microbenchmark(
  
  meth1 = {
    beta <- as.matrix(crossprod(X, y))
    M <- y
    M@x = M@x - suvC(X, t(beta), irow, pcol)
    M <- naive_MC(as.matrix(M))
  },
  meth2 = {
    Xterms = utils$GetXterms(X)
    Y_naive = naive_MC(as.matrix(y))
    beta = Xterms$X1 %*% Y_naive
    Xbeta <- X %*% beta   #svdH$u %*% (svdH$v  %*% Y_naive)
    M <- Y_naive - Xbeta
  },
  times = 100L
  
)