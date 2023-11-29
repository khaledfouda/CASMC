
# Y_data : output matrix. a matrix without missing data
# X_data : row features data
# Z_data : col features data
# mask :  mask matrix
# Y_pred_x : predicted output using Mao's method on the row covariates
# Y_pred_z : predicted output using Mao's method on the col covariates
# Y_pred_new <- Y_pred_x + Y_pred_z : predicted output

# for functions "SMCfit_cv", and "SMCfit", see the file "SMC_functions6.R" or Mao's github.

#alpha <- 1  # the agreement penalty parameter
#max_iter <- 5
print(paste0("Agreement Penalty = ", alpha))
set.seed(2023)
start_time = Sys.time()
#--- The proposed method alternates between the application of Mao's method on the row covariates and on the column covariates. 
#--- Hence, first we try to decide whether to start with row covairates or column covariates.
error_function <- function(preds, orig){
  return(sum( (preds-orig)^2 )/ sum((orig)^2))
}

##---- We run Mao's method on the row covariates
cv_ratio_SMC_x = SMCfit_cv(Y_data, X_data, mask, nfolds = 5, 
                           tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                           alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
fit_mao_x <- SMCfit(Y_data*mask, X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)
Y_pred_x <- fit_mao_x$Ahat
##--- We run Mao's method on the column covariates
###--- Attention : Pay attention to the dimension of the quantities, 
###--- use the transpose operator when switching from the row covariates to the column covariates.
t_Y_data <- t(Y_data)
t_mask <- t(mask)
cv_ratio_SMC_z = SMCfit_cv(t_Y_data, Z_data, t_mask, nfolds = 5, 
                           tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                           alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
fit_mao_z <- SMCfit(t_Y_data*t_mask, Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)
Y_pred_z <- t(fit_mao_z$Ahat)

test_indices <- which(W_data==0, arr.ind = TRUE)
error_x <- error_function(Y_pred_x[test_indices],Y_data[test_indices])
#sqrt(mean((Y_data - Y_pred_x)[test_indices]**2))
error_z <- error_function(Y_pred_z[test_indices],Y_data[test_indices])
  #sqrt(mean((Y_data - Y_pred_z)[test_indices]**2))
Y_combined = (Y_pred_x + Y_pred_z)/2
error_combined <- error_function(Y_combined[test_indices],Y_data[test_indices])
  #sqrt(mean((Y_data - Y_combined)[test_indices]**2))
print(paste0("MSE X=",error_x,", MSE Z=",error_z, ", MSE Avg=",error_z))

##---  If using first the row covariates (x) minimizes the cv-error then we start with the row covariates 
if(cv_ratio_SMC_x$cv_min < cv_ratio_SMC_z$cv_min){
  
  ###--- first we run one iteration to initialize all the quantities
  output_z <- t_Y_data - (1-alpha)*t(Y_pred_x)
  cv_ratio_SMC_z = SMCfit_cv(output_z/(1+alpha), Z_data, t_mask, nfolds = 5, 
                             tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                             alpha_grid = seq(0.992,1, length = 10))
  fit_mao_z <- SMCfit(output_z*mask/(1+alpha), Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)
  Y_pred_z <- t(fit_mao_z$Ahat)
  
  output_x <- Y_data - (1-alpha)*Y_pred_z
  cv_ratio_SMC_x = SMCfit_cv(output_x/(1+alpha), X_data, mask, nfolds = 5, 
                             tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                             alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
  fit_mao_x <- SMCfit(output_x*mask/(1+alpha), X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)
  Y_pred_x <- fit_mao_x$Ahat
  Y_pred_new <- Y_pred_x + Y_pred_z
  
  ###--- we iterate 100 times or until convergence 
  for(iter in 1:max_iter){
    Y_pred_old <- Y_pred_new
    output_z <- t_Y_data - (1-alpha)*t(Y_pred_x)
    cv_ratio_SMC_z = SMCfit_cv(output_z/(1+alpha), Z_data, t_mask, nfolds = 5, 
                               tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                               alpha_grid = seq(0.992,1, length = 10))
    fit_mao_z <- SMCfit(output_z*mask/(1+alpha), Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)
    Y_pred_z <- t(fit_mao_z$Ahat)
    
    output_x <- Y_data - (1-alpha)*Y_pred_z
    cv_ratio_SMC_x = SMCfit_cv(output_x/(1+alpha), X_data, mask, nfolds = 5, 
                               tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                               alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
    fit_mao_x <- SMCfit(output_x*mask/(1+alpha), X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)
    Y_pred_x <- fit_mao_x$Ahat
    Y_pred_new <- Y_pred_x + Y_pred_z
    diff <- sqrt(mean((Y_pred_new-Y_pred_old)**2))
    test_indices <- which(W_data==0, arr.ind = TRUE)
    mao_error_coop <- Y_data - Y_pred_new
    mao_error_coop <- mao_error_coop[test_indices]
    mao_error_coop <- sqrt(mean(mao_error_coop**2))
    mao_error_coop
    print(paste0("iter=",iter,", diff=",diff, ", error=",mao_error_coop))
    if(diff< 1e-5){break}
    
  }
  
}else{
  ##---  If using first the col covariates (z) minimizes the cv-error then we start with the col covariates  
  ###--- first we run one iteration to initialize all the quantities
  output_x <- Y_data - (1-alpha)*Y_pred_z
  cv_ratio_SMC_x = SMCfit_cv(output_x/(1+alpha), X_data, mask, nfolds = 5, 
                             tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                             alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
  fit_mao_x <- SMCfit(output_x*mask/(1+alpha), X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)
  Y_pred_x <- fit_mao_x$Ahat
  
  output_z <- t_Y_data - (1-alpha)*t(Y_pred_x)
  cv_ratio_SMC_z = SMCfit_cv(output_z/(1+alpha), Z_data, t_mask, nfolds = 5, 
                             tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                             alpha_grid = seq(0.992,1, length = 10))
  fit_mao_z <- SMCfit(output_z*mask/(1+alpha), Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)
  Y_pred_z <- t(fit_mao_z$Ahat)
  Y_pred_new <- Y_pred_x + Y_pred_z
  
  ###--- we iterate 100 times or until convergence 
  for(iter in 1:max_iter){
    Y_pred_old <- Y_pred_new
    output_x <- Y_data - (1-alpha)*Y_pred_z
    cv_ratio_SMC_x = SMCfit_cv(output_x/(1+alpha), X_data, mask, nfolds = 5, 
                               tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                               alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.
    fit_mao_x <- SMCfit(output_x*mask/(1+alpha), X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)
    Y_pred_x <- fit_mao_x$Ahat
    
    output_z <- t_Y_data - (1-alpha)*t(Y_pred_x)
    cv_ratio_SMC_z = SMCfit_cv(output_z/(1+alpha), Z_data, t_mask, nfolds = 5, 
                               tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                               alpha_grid = seq(0.992,1, length = 10))
    fit_mao_z <- SMCfit(output_z*mask/(1+alpha), Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)
    Y_pred_z <- t(fit_mao_z$Ahat)
    Y_pred_new <- Y_pred_x + Y_pred_z
    diff <- sqrt(mean((Y_pred_new-Y_pred_old)**2))
    test_indices <- which(W_data==0, arr.ind = TRUE)
    #mao_error_coop <- Y_data - Y_pred_new
    #mao_error_coop <- mao_error_coop[test_indices]
    #mao_error_coop <- sqrt(mean(mao_error_coop**2))
    mao_error_coop <- error_function(Y_pred_new[test_indices],Y_data[test_indices])
    print(paste0("iter=",iter,", diff=",diff, ", error=",mao_error_coop))
    if(diff< 1e-5){break}
  }
  
}


time_taken <- round(as.numeric(difftime(Sys.time(), start_time, units='mins')),3)   

print(paste0("Run time is ",time_taken, " minutes."))
print(paste0("Average Run time per iteration is ",time_taken/iter, " minutes."))

#--- code to get the RMSE
test_indices <- which(W_data==0, arr.ind = TRUE)
mao_error_coop <- Y_data - Y_pred_old
mao_error_coop <- mao_error_coop[test_indices]
mao_error_coop <- sqrt(mean(mao_error_coop**2))
mao_error_coop




