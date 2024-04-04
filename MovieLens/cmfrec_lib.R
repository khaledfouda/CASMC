library(cmfrec)
library(Matrix)
library(MatrixExtra)
setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP")
suppressMessages(suppressWarnings(source("./code_files/import_lib.R", local = FALSE)))
source("./MovieLens/load_data.R")
#---------------

dat <- load_movielens_100k("a")

X_train <- as.coo.matrix(dat$train_unscaled)
str(X_train)
X_test <- as.coo.matrix(dat$test)
str(X_test)

model.classic <- CMF(X_train, k=25, lambda=0.1, scale_lam=TRUE, verbose=FALSE)
