setwd("/mnt/campus/math/research/kfouda/main/HEC/Youssef/HEC_MAO_COOP/code_files/AIS")

source("./AIS_fit.R")
source("./Approximate_SVT.R")
source("../Mao_import_lib.R")



gen.dat <- generate_simulation_data_ysf(2,600,600,10,10, missing_prob = 0.9,coll=TRUE)

# W_valid <- matrix.split.train.test(gen.dat$W, testp=0.2)
# Y_train = (gen.dat$Y * W_valid)
# Y_valid = gen.dat$Y[W_valid==0]


fit_ = AIS_fit(gen.dat$Y, lambda=30, trace.it = TRUE)
Y = gen.dat$Y; lambda=30; maxR = 10; maxit=100; tol=1e-6; thresh=1e-3;decay=0.8; trace.it=FALSE
