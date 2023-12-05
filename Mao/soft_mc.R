

S_data <- Y_data
S_data[W_data==0]=NA 

###uses regular matrix method for matrices with NAs


fit_soft=softImpute(S_data, rank=50, lambda=30)

test_indices <- which(W_data==0, arr.ind = TRUE)
soft_error <- complete(S_data, fit_soft)
soft_error <- Y_data - soft_error
soft_error <- soft_error[test_indices]
# soft_error <- impute(fit_soft, i=test_indices[,1], j=test_indices[,2])
soft_error <- sqrt(mean(soft_error**2))
# soft_error
