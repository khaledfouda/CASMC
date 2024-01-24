

cv_ratio_SMC_x = SMCfit_cv(Y_data, X_data, mask, nfolds = 5, 
                         tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                         alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.

fit_mao_x <- SMCfit(Y_data*mask, X_data, cv_ratio_SMC_x$tau_beta_ratio, cv_ratio_SMC_x$tau_svd_ratio, cv_ratio_SMC_x$alpha_ratio)

Y_pred_x <- fit_mao_x$Ahat

test_indices <- which(W_data==0, arr.ind = TRUE)
mao_error_x <- Y_data - Y_pred_x
mao_error_x <- mao_error_x[test_indices]
mao_error_x <- sqrt(mean(mao_error_x**2))
mao_error_x 
