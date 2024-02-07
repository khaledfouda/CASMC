
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
X <- gen.dat$X
lambda2 = 31
max.rank = 15
start_time <- Sys.time()
beta_partial = solve(t(gen.dat$X) %*% gen.dat$X) %*% t(gen.dat$X)
set.seed(2023);sout <- simpute.als.cov(Y_train, gen.dat$X, beta_partial,J = max.rank, thresh =  1e-6,
                        lambda= lambda2,trace.it = T,warm.start = NULL, maxit=100)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))
sout$A_hat = sout$u %*% (sout$d * t(sout$v))
fiti <- sout
v=as.matrix(fiti$v)
vd=v*outer(rep(1,nrow(v)),fiti$d)

fiti=fits2

sout$A_hat = fiti$u %*% t(vd)  + X %*% fiti$beta.estim

print(paste("Test error =", round(test_error(sout$A_hat[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (sout$A_hat[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))
print(paste("Test error =", round(test_error(sout$u %*% (sout$d * t(sout$v)), gen.dat$B),5)))

xbeta.orig = X %*% fiti$beta.estim
#-----------------------
# new model
y = Y_train
thresh=1e-6
y[y==0] = NA
ys <- as(y, "Incomplete")
yvalid = gen.dat$Y
yvalid[W_valid==1] = NA
yvalid[yvalid==0] = NA
yvalid <- as(yvalid, "Incomplete")
#$$
#y=ys; X=gen.dat$X; H=NULL; J = 2; thresh = 1e-05; lambda=2; 
#maxit=100;trace.it=T;warm.start=NULL;final.svd=FALSE; patience=3
#$$

#--H
H = X %*% solve(t(X) %*% X) %*% t(X)
svdH <- fast.svd(H, thresh)
J_H <- sum(svdH$d > 1e-6)
print(J_H)
svdH$u = svdH$u[,1:J_H]
svdH$v = svdH$d[1:J_H] * t(svdH$v[,1:J_H])
svdH$d = NULL

#------------------------
start_time <- Sys.time()
set.seed(2023);fits <- simpute.als.fit_splr(y=ys, yvalid=yvalid, X=gen.dat$X,  trace=F, J=31,
                                                    thresh=1e-6, lambda=31, H=NULL,svdH=svdH,
                                   final.svd = T,maxit = 300, patience=1)
print(paste("Execution time is",round(as.numeric(difftime(Sys.time(), start_time,units = "secs")),2), "seconds"))


M = fits$u %*% (fits$d * t(fits$v))
xbeta.S = ys; xbeta.S@x = fits$xbeta.obs
M.S = M; M.S[Y_train==0] = 0

print(paste("Test error =", round(test_error(M, gen.dat$B),6)))

sum(round(xbeta.orig[Y_train!=0],2) == round(xbeta.S[xbeta.S!=0],2))/sum(Y_train!=0)

preds = H %*% (Y_train - M.S) + M

sqrt(mean( (M[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))

preds =  H%*% (-M) +M
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))

preds = (xbeta.S + M)[ys!=0]
print(paste("Test error =", round(test_error(preds, gen.dat$A[Y_train!=0]),5)))
sqrt(mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))

set.seed(2023);fits2 <- simpute.als.cov(Y_train, gen.dat$X, beta_partial,J = max.rank, thresh =  1e-6,
                                       lambda= lambda2,trace.it = T,warm.start = fits.out, maxit=100)

fits.out = list(u=fits$u, d=fits$d, v=fits$v, beta.estim=fits$beta_estim)
ytmp = Y_train 
ytmp[Y_train==0] = (M)[Y_train==0]
fits$beta_estim = beta_partial %*% (ytmp)
fits$beta_estim = H %*% (Y_train-M)
dim(y)
dim(fits$beta_estim)

preds <- predsa <- M + fits$beta_estim

best_score = Inf
patience = 3
counter = 0
alpha = 0.1
H_part = (Diagonal(800) + alpha *  H )
for(i in 1:400){
   if(i == 1) preds <- predsa <-  fits$beta_estim
preds <- H_part %*% (preds)
RMSE = round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)
#mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 )
if(RMSE < best_score){
   best_score = RMSE
   best_i = i
   best_preds = preds
   counter = 0
   #print(paste(i,RMSE))
}else counter = counter + 1
if(counter > patience){
   print(paste("Exit on",i))
   break
}
}
print(paste(best_i, best_score))
preds = best_preds + M
print(paste("Test error =", round(test_error(preds[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))
sqrt(mean( (preds[gen.dat$W==0]-gen.dat$A[gen.dat$W==0])^2 ))

ytmp[is.na(y)] = (M - H %*% M)[is.na(y)]

fits$beta_estims = H %*% (gen.dat$Y -  fits$u %*% (fits$d * t(fits$v)))
all(round(fits$beta_estims,10) == 0)
preds = M+ fits$beta_estim #fits$beta_estim + M
preds = M +  H %*% fits$M_sum #best_preds + M
preds = M + (H %*% ( 1 *ys + 10 *M))
preds = M + H %*% (ys - M)




preds = M + best_preds
dim(fiti$beta.estim)
dim(X)
#----------
fits <- NULL
start_time <- Sys.time()
lambda2_vals = seq(140,0,-5)
scores = rep(NA, length(lambda2_vals))
for(i in 1:length(lambda2_vals)){
   
set.seed(2023);fits <- simpute.als.fit_splr(ys, yvalid, gen.dat$X, trace=F, J=5,
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





#----------------------
xb =fits$xbeta.obs


length(xb)

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


soft_estim = complete(Y_train, fits)
print(paste("Test error =", round(test_error(soft_estim[gen.dat$W==0], gen.dat$A[gen.dat$W==0]),5)))

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
preds <- H %*% (ys-preds) + preds

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
###########################################################

n = 800
m = 900


H = matrix(rnorm(n*n), n,n)
mask = Matrix( rbinom(n*m,1,0.1),n,m )
S = matrix(rnorm(n*m), n,m) * mask
S = Matrix(S, sparse = TRUE)
mask = matrix(1, n, n)   

sp = S@p
si = S@i
sx = S@x
r = length(sx)

sparse_prod_R = function(H, sp, si, sx, n, m, r){


index = 1
result = numeric(r)

for(j in 1:m){
   
   jstart = sp[j] +1
   jend = sp[j+1]
   
   if(jstart > jend) next
   
   ind = jstart:jend
   # observed subset for column j in S
   si_sub = si[ind] + 1
   
   # we now compute the subset of rows in H which correspond to observed rows in S
   
   sx_sub = sx[ind]
   # to be multiplied by each column in H
   
   h_sub = H[,si_sub]
   
   
   for(hrow in si_sub){
      result[index] = geometry::dot(sx_sub, h_sub[hrow,])
      index = index + 1
   }
}
return(result)
}


system.time({
   for(i in 1:30) 
         result = sparse_prod_R(H, sp, si, sx, n, m, r)
      
      })

length(sx)
sp[1:10]

system.time({
   for(i in 1:30) 
      
result2 = sparse_prod(n, m, r, H, sp, si, sx)
})


system.time({
   for(i in 1:30) 
      results0 = (svd_result$u %*% ((svd_result$d * t(svd_result$v)) %*% S) )[ss]
   #result0 = (H %*% S)[ss]
})









sparse_multiply <- function(H, Y) {
   result <- matrix(0, nrow = nrow(H), ncol = ncol(Y))
   for (i in 1:nrow(Y)) {
      non_zero_indices <- which(Y[i, ] != 0)
      if (length(non_zero_indices) > 0) {
         result[i, non_zero_indices] <- H[i, ] %*% Y[i, non_zero_indices]
      }
   }
   return(result)
}

# Perform the sparse multiplication
result <- sparse_multiply(H, S)







storage.mode(H) = "double"
storage.mode(sx) = "double"
storage.mode(si) = "integer"
storage.mode(sp) = "integer"
storage.mode(n) = "integer"
storage.mode(m) = "integer"
storage.mode(r) = "integer"

# Preallocate the result vector
result = double(r)
system.time(
   {
   for(i in 1:30)
# Call the Fortran subroutine
result2 = .Fortran("sparse_prod",
         as.integer(n), as.integer(m), as.integer(r),
         as.double(H), as.integer(sp), as.integer(si),
         as.double(sx), result = as.double(result)
)$result

   }
)


ss = S!=0
S2 <- as(S, "Incomplete")


result0[1:4]
result[1:4]
result2[1:4]


all(round(result0,5) == round(result2,5))

h_sub[hrow,]
sx_sub

all(round((H %*% S)[S!=0],5) == round(result,5))

length((H %*% S)[S!=0])
result %>% length
result[1:10]
(H %*% S)[S!=0][1:10]

length(S@x)

library(geometry)
install.packages("geometry")


k <- 10  # Choose a suitable rank
svd_result <- irlba(H, nv = k)

# Reconstruct lower-rank approximation of H
H_approx <- svd_result$u %*% (svd_result$d * t(svd_result$v))


dim(H_approx)
