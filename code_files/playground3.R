
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")


gen.dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.9,coll=T)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]
#X <- gen.dat$X
X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
beta_partial = MASS::ginv(t(X_r$X) %*% X_r$X) %*% t(X_r$X)
#--
y = yfill = Y_train
y[y==0] = NA
ys <- as(y, "Incomplete")
xbeta.sparse = ys
# yvalid = gen.dat$Y
# yvalid[W_valid==1] = NA
# yvalid[yvalid==0] = NA
# yvalid <- as(yvalid, "Incomplete")
#---
lambda2 = 11.08761
max.rank = 5
#---

start_time <- Sys.time()
set.seed(2023);fits <- simpute.als.fit_splr(y=ys, svdH=X_r$svdH,  trace=F, J=max.rank,
                                            thresh=1e-6, lambda=lambda2, return_obj = F, init = "naive",
                                            final.svd = T,maxit = 100, warm.start = NULL,Px = beta_partial)

print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fits$M = fits$u %*% (fits$d * t(fits$v))

xbeta =  X_r$X %*% fits$beta.obs
xbeta.orig = gen.dat$X %*% gen.dat$beta.x

fits$xbeta.obs[1:8] %>% round(4)
xbeta[Y_train==0][1:8]  %>% round(4)
init_xbeta[Y_train==0][1:8]%>% round(4)
xbeta.orig[Y_train==0][1:8]%>% round(4)
new.xbeta[Y_train==0][1:8]%>% round(4)
(gen.dat$X %*% soutl$beta_hat)[Y_train==0][1:8]%>% round(4)

fits$beta.obs[1:5,1:5]
soutl$beta_hat[1:5,1:5]

test_error(fits$beta.obs, gen.dat$beta.x)
test_error(X_r$X %*%fits$beta.obs, gen.dat$X %*%soutl$beta_hat)
test_error(fits$xbeta.obs, (gen.dat$X %*% soutl$beta_hat)[Y_train!=0] ) # better than new.xbeta
test_error(new.xbeta[Y_train!=0], (gen.dat$X %*% soutl$beta_hat)[Y_train!=0] )



test_error(fits$xbeta.obs, (gen.dat$X %*% gen.dat$beta.x)[Y_train!=0] ) # better than new.xbeta and soutl
test_error(new.xbeta[Y_train==0], (gen.dat$X %*% gen.dat$beta.x)[Y_train==0] )
test_error((gen.dat$X %*% soutl$beta_hat)[Y_train==0], (gen.dat$X %*% gen.dat$beta.x)[Y_train==0] )
test_error(init_xbeta[Y_train==0], (gen.dat$X %*% gen.dat$beta.x)[Y_train==0] )


test_error(soutl$beta_hat, gen.dat$beta.x)
test_error(ginv(X_r$X) %*% init_xbeta, gen.dat$beta.x)
test_error(ginv(X_r$X) %*% new.xbeta, gen.dat$beta.x)
test_error(ginv(X_r$X) %*% new.xbeta, soutl$beta_hat)

fits$A = fits$M + X_r$X %*% fits$beta.obs
#---
#xbeta.sparse = as.matrix(xbeta.sparse)
#xbeta.sparse[is.na(xbeta.sparse)] = 0

xbeta.sparse@x = fits$xbeta.obs
init_xbeta = naive_MC(as.matrix(xbeta.sparse))
#yfill[Y_train==0] = fits$M[Y_train==0]
#init_xbeta =  X_r$svdH$u %*% (X_r$svdH$v %*% yfill)
#init_xbeta =  as.matrix(X_r$X %*% fits$beta.obs)#X_r$svdH$u %*% (X_r$svdH$v %*% yfill)
warm.start = propack.svd(init_xbeta, X_r$rank)
start_time <- Sys.time()
fits3 = simpute.als.splr.fit.nocov.fixedJ(xbeta.sparse, X_r$rank, maxit=200, final.trim = F,thresh=1e-5,
                                          warm.start = warm.start, trace.it=T, return_obj = F)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

new.xbeta = fits3$u %*% (fits3$d * t(fits3$v))
fits$A = fits$M + fits3$u %*% (fits3$d * t(fits3$v))
fits$beta_hat = MASS::ginv(X_r$X) %*% fits3$u %*% (fits3$d * t(fits3$v))
print(paste("Test error =", round(test_error(fits$beta_hat, gen.dat$beta.x),5)))

sqrt(mean( ((fits$A)[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
print(paste("Test error =", round(test_error((fits$A)[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
print(paste("Test error =", round(test_error(fits$M, gen.dat$B),5)))
#----------
# New fit function
start_time <- Sys.time()
fit4 = simpute.als.splr.fit.beta(as.matrix(xbeta.sparse), X_r$X, ginv(X_r$X), X_r$J,
                                 trace.it = F,final.trim = F)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
beta_hat4 = ginv(X_r$X) %*% fit4$K
print(paste("Test error =", round(test_error(beta_hat4, gen.dat$beta.x),5)))


#---------
#fits2 = restore.beta(Y_train, fits$M, xbeta.sparse, X_r$rank, thresh=1e-6, maxit=300, trace.it=T)
sum(round(sout$M,3)==round(fits$M,3)) / length(sout$M)


all(xbeta.sparse[!is.na(xbeta.sparse)] == fits$xbeta.obs)
sum(round(sout$Xbeta[Y_train!=0],2) == round(fits$xbeta.obs,2)) / length(fits$xbeta.obs)

dim(MASS::ginv(gen.dat$X))

sum(round(sout$beta.estim,2) == round(beta,2)) / length(beta)
fits3$J
###################################################
# cross-validation
set.seed(2023)
gen.dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.8,coll=F)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]
X_r = reduced_hat_decomp(gen.dat$X, 1e-2)
start_time <- Sys.time()
fitcv = simpute.cov.cv_splr(Y_train, X_r, Y_valid, W_valid, trace=T)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

fitcv$A = fitcv$B_hat + fitcv$xbeta$u %*% (fitcv$xbeta$d * t(fitcv$xbeta$v))
fitcv$beta = (MASS::ginv(X_r$X) %*% fitcv$xbeta$u) %*% (fitcv$xbeta$d * t(fitcv$xbeta$v))
sqrt(mean( ((fitcv$A)[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
print(paste("Test error =", round(test_error((fitcv$A)[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
print(paste("Test error =", round(test_error(fitcv$B_hat, gen.dat$B),5)))
print(paste("Test error =", round(test_error(fitcv$beta, gen.dat$beta.x),5)))
fitcv$lambda
qr(fitcv$A)$rank
#------
####################################
# K- Fold
start_time <- Sys.time()
fitkf = simpute.cov.Kf_splr(gen.dat$Y, X_r, gen.dat$W, beta_partial, 3, trace=T)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))

print(paste("Test error =", round(test_error((fitkf$A_hat)[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
print(paste("Test error =", round(test_error(fitkf$B_hat, gen.dat$B),5)))
print(paste("Test error =", round(test_error(MASS::ginv(X_r$X)%*%fitkf$beta_hat, gen.dat$beta.x),5)))
fitkf$lambda2
fitkf$rank_A

###################################################
#' Best case scenario
#' K-fold with ALS
sout <- simpute.cov.kfold(gen.dat$Y, gen.dat$X, gen.dat$W, n_folds = 3, print.best = FALSE,
                          trace=TRUE, rank.limit = 30, lambda1=0,n1n2 = 1, warm=NULL,tol = 2)
soutl <- simpute.cov.kfold.lambda1(gen.dat$Y, gen.dat$X, gen.dat$W, sout$lambda2, n_folds = 3, print.best = FALSE, 
                                  trace=TRUE,lambda1.grid = seq(0,20,length.out=20) ,n1n2 = 1, warm=NULL,
                                  J=c(sout$J))


print(paste("Test error =", round(test_error((soutl$A_hat)[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
print(paste("Test error =", round(test_error(soutl$B_hat, gen.dat$B),5)))
print(paste("Test error =", round(test_error(soutl$beta_hat, gen.dat$beta.x),5)))
soutl$lambda1
soutl$lambda2



######################################################
# analyze obj
# plot(1:200, fits3$obj, pch=4, col="blue")

fits_obj = c()
for(i in 1:15){
   fitss = simpute.als.splr.fit.nocov.fixedJ(xbeta.sparse, X_r$rank, maxit=200, final.trim = T,
                                             warm.start = NULL, trace.it=F, return_obj = T)
   fits_obj = c(fits_obj, fitss$obj)
   #points(1:200, fitss$obj, pch=".")
}

data <- data.frame(iterations = rep(1:200, 15),
                   objective = fits_obj,
                   fit_i = rep(1:15, each=200))





lowest_points <- data %>%
   group_by(fit_i) %>%
   summarise(iterations = which.min(objective), min_objective = min(objective,na.rm = T)) %>%
   ungroup()

ggplot(data, aes(x = iterations, y = objective, color = fit_i)) +
   geom_line() +  
   geom_point(data = lowest_points, aes(x = iterations, y = min_objective, color = fit_i),
              size = 3, shape = 23, fill = "white") +

   geom_line(data=data.frame(x=1:200,y=fits3$obj, fit_i=0),aes(x,y)) +
   #  geom_text(data = lowest_points, 
#             aes(x = iterations, y = min_objective, label = paste("Min:", round(min_objective, 2))),
#             nudge_y = -5, color = "black", size = 3) +
   theme_minimal() +
   labs(title = "Convergence of Fitting Methods",
        subtitle = "Each line represents a fitting method's convergence over iterations.\nThe lowest points are highlighted.",
        x = "Iteration",
        y = "Objective Value",
        color = "Method") +
   theme(legend.position = "bottom",
         plot.title = element_text(hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5),
         legend.title.align = 0.5)


fits3$obj[1:fits3$iter]


points(1:300, fits3$obj, pch=".")

#----------------
start_time <- Sys.time()

set.seed(2023);sout <- simpute.als.cov(Y_train, X_r$X, beta_partial,J = max.rank, thresh =  1e-6,
                                       lambda= lambda2,trace.it = T,warm.start = NULL, maxit=100)

print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
sout$M = sout$u %*% (sout$d * t(sout$v))
sout$Xbeta = X %*% sout$beta.estim
sout$A = sout$M  + sout$Xbeta

print(paste("Test error =", round(test_error(sout$A[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (sout$A[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
print(paste("Test error =", round(test_error(sout$M, gen.dat$B),5)))
print(paste("Test error =", round(test_error(sout$beta.estim, gen.dat$beta.x),5)))




#---------------------------------------------------------------------------------------------------------


# Cross-validation

best_fit = simpute.cov.cv_splr_no_patience(Y_train, X_reduced$svdH, Y_valid, W_valid,warm = best_fit$best_fit,
                                            trace = F, rank.limit=30,rank.step=4,patience = 1,
                                           rank.init = 2, lambda.factor = 1/2, n.lambda = 30)
best_fit = simpute.cov.Kf_splr_no_patience_v2(gen.dat$Y, X_reduced$svdH, gen.dat$W, n_folds=10,
                                            trace = F, rank.init=2, lambda.factor = 1/4,n.lambda = 30,
                                            rank.limit=30,rank.step=4,patience = 1)


yfill = gen.dat$Y
fits = best_fit$best_fit
best_fit$B_hat = fits$u %*% (fits$d * t(fits$v))
yfill[gen.dat$Y==0] = (best_fit$B_hat)[gen.dat$Y==0]
fits.out = list(u=fits$u, d=fits$d, v=fits$v, beta.estim=beta_partial %*% yfill)
beta_partial = MASS::ginv(t(X_reduced$X) %*% X_reduced$X) %*% t(X_reduced$X)


set.seed(2023);fits2 <- simpute.als.cov(gen.dat$Y, X_reduced$X, beta_partial,J = best_fit$rank_B, thresh =  1e-6,
                                        lambda= best_fit$lambda,trace.it = F,warm.start = fits.out, maxit=100)

fits2$M = fits2$u %*% (fits2$d * t(fits2$v))
fits2$A_hat = fits2$M  + X_reduced$X %*% fits2$beta.estim
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
