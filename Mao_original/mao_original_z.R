

t_Y_data <- t(Y_data)
t_mask <- t(mask)
cv_ratio_SMC_z = SMCfit_cv(t_Y_data, Z_data, t_mask, nfolds = 5, 
                         tau1_grid = seq(0, 2, length = 20), tau2_grid = seq(0.9, 0.1, length = 20),
                         alpha_grid = seq(0.992,1, length = 10))  # Choose the tuning parameter by cross-validation.

fit_mao_z <- SMCfit(t_Y_data*t_mask, Z_data, cv_ratio_SMC_z$tau_beta_ratio, cv_ratio_SMC_z$tau_svd_ratio, cv_ratio_SMC_z$alpha_ratio)

Y_pred_z <- fit_mao_z$Ahat

test_indices <- which(t_mask==0, arr.ind = TRUE)
mao_error_z <- t_Y_data - Y_pred_z
mao_error_z <- mao_error_z[test_indices]
mao_error_z <- sqrt(mean(mao_error_z**2))
mao_error_z 
