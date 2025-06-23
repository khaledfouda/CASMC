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

VDsq = t(Dsq * t(V))

all(
  t(crossprod(y, U) + VDsq) ==
    B
)

res1  <- microbenchmark(
  old_D = {Dsq / (Dsq + lambda.M)},
  new_D = {1/(1 + lambda.M * (1/Dsq) )},
  times = 1000L, units = "ms"
)
print(res1)
plot(res1)


Dstar = (Dsq / (Dsq + lambda.M))

B = t(U) %*% y + t(VDsq)
B = as.matrix(t((B) * Dstar))


B2 = (crossprod(y, U) + VDsq)
B2 = as.matrix(B2 %*% diag(Dstar))
B2 = sweep(B2, 2L, (Dsq / (Dsq + lambda.M)), `*`)

B2[] <- B2

all(B==(B2))

B

res2  <- microbenchmark(
  old_D = {
    #B = t(U) %*% y + t(VDsq)
    B22 = (B2 %*% diag(Dstar))
    B22 = as.matrix(B22)
    #B11 = as.matrix(t((B) * Dstar))
  },
  new_D = {
    #B2 = (crossprod(y, U) + VDsq)
    B22 = as.matrix(B2 %*% diag(Dstar))
    # B22 = sweep(B2, 2L, Dstar, `*`)
  },
  times = 1000L#, units = "sec"
)
print(res2)
plot(res2)
fast.svd(B22)
t(B2)