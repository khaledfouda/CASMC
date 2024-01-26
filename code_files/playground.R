
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
source("./code_files/import_lib.R")

gen.dat <- generate_simulation_data_ysf(1,500,500,5,10, missing_prob = 0.9,coll=T,seed=3023)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_valid <- gen.dat$Y[W_valid==0]

Y_train <- gen.dat$Y
Y_train[Y_train==0] = NA
y = Y_train %>% as.matrix()
ys <- as(y, "Incomplete")

ysn <- y



onesSparse <- ys
onesSparse@x[] <- 1
ysub <- as(y, "Incomplete")
resultSparse <- ysub * onesSparse
resultSparse@
resultMat <- as.matrix(resultSparse)

fits <- simpute.als.fit_Incomplete(ys, gen.dat$X, trace=TRUE, J=8, lambda=3,final.svd = F)

preds <- fits$u %*% (fits$d * t(fits$v))
test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0])

fits$rank

y=ys; X=gen.dat$X; J = 2; thresh = 1e-05; lambda=2; 
maxit=100;trace.it=T;warm.start=NULL;final.svd=FALSE


fits$rank
fits$lambda
system.time({fitss <- simpute.als.fit_Incomplete(ys, J = 3, lambda=1.9)})


#---------------------------------
# old model 
gen.dat <- generate_simulation_data_ysf(2,800,800,10,10, missing_prob = 0.9,coll=FALSE)
W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
Y_train = (gen.dat$Y * W_valid)
Y_valid = gen.dat$Y[W_valid==0]

lambda2 = 31
max.rank = 15
start_time <- Sys.time()
beta_partial = solve(t(gen.dat$X) %*% gen.dat$X) %*% t(gen.dat$X)
set.seed(2023);sout <- simpute.als.cov(Y_train, gen.dat$X, beta_partial,J = max.rank, thresh =  1e-6,
                        lambda= lambda2,trace.it = T,warm.start = NULL, maxit=100)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
sout$A_hat = sout$u %*% (sout$d * t(sout$v))
print(paste("Test error =", round(test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (sout$A_hat[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))

#-----------------------
# new model
y = Y_train
y[y==0] = NA
ys <- as(y, "Incomplete")
yvalid = gen.dat$Y
yvalid[W_valid==1] = NA
yvalid[yvalid==0] = NA
yvalid <- as(yvalid, "Incomplete")

start_time <- Sys.time()
set.seed(2023);fits <- simpute.als.fit_Incomplete_2(ys, gen.dat$X, yvalid, trace=T, J=max.rank,
                                                    thresh=1e-6, lambda=31,
                                   final.svd = T,maxit = 300, patience=1)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
preds <- fits$u %*% (fits$d * t(fits$v))
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
#----------
fits <- NULL
start_time <- Sys.time()
lambda2_vals = seq(140,0,-5)
scores = rep(NA, length(lambda2_vals))
for(i in 1:length(lambda2_vals)){
   
set.seed(2023);fits <- simpute.als.fit_Incomplete_2(ys, gen.dat$X, yvalid, trace=F, J=5,
                                                    thresh=1e-6, lambda=lambda2_vals[i], warm.start = fits,
                                                    final.svd = T,maxit = 1000, patience=1)
preds <- fits$u %*% (fits$d * t(fits$v))
scores[i] = round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)#fits$best_score
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))

}

min(scores)
lambda2_vals[which.min(scores)]

print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
sqrt(mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))







#-------------------------------
timespent = rep(0,2)
start_time <- Sys.time()
for(i in 1:1000){
#hat2 <- X %*% solve(t(X) %*% X) %*% t(X)}
 p = ncol(X)
 Q <- qr.Q(Matrix::qr(X))[,1:p]
 hat <- Q %*% t(Q)}
#svdX <- fast.svd(X)
#hat = svdX$u %*% Diagonal(length(svdX$d)) %*% t(svdX$u)
# onesSparse <- ys
# onesSparse@x[] <- 1
# hat_sp <- as(as.matrix(hat), "Incomplete")
# hat_sp <- hat_sp * onesSparse
timespent[1] = timespent[1] + as.numeric(difftime(Sys.time(), start_time,units = "secs"))
timespent
#---------------------------
# new model
y = Y_train
y[y==0] = NA
ys <- as(y, "Incomplete")
start_time <- Sys.time()
fits <- simpute.als.fit_Incomplete_old(ys, trace=F, J=max.rank, thresh=1e-3, lambda=lambda2,
                                   final.svd = T)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
preds <- fits$u %*% (fits$d * t(fits$v))
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))

#---------------------------------------------
# new model
y = Y_train
y[y==0] = NA
ys <- as(y, "Incomplete")
start_time <- Sys.time()
fits <- softImpute(ys, trace=F, rank.max=max.rank, thresh=1e-3, lambda=lambda2,
                                       final.svd = T)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
preds <- fits$u %*% (fits$d * t(fits$v))
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))

#-------------------------------------------------------------------------------------------------
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
#ysc <- biScale(ys, col.scale=FALSE,row.scale=FALSE,trace=TRUE)
system.time({fitss <- softImpute(ys, rank.max = 3, lambda=1.9)})
system.time({fits <- softImpute(Y_train, rank.max = 3, lambda=1.9)})

system.time({complete(ys, fitss)})
# fitss <- deBias(ysc, fitss)
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

ys -> mat
rows <- mat@i + 1
cols <- rep(1:ncol(mat), diff(mat@p))
ynew <- y
# Combine into a matrix where each row is a (row, column) pair
indices <- cbind(rows, cols)

dim(indices)
sum(!is.na(y))

y[indices[,1],indices[,2]]

mask <- ys != 0
ynew[mask] 

sum(ys!=0)

# Optional: Convert to a list of indices
indices_list <- split(indices, seq(nrow(indices)))
