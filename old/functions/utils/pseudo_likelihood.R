utils$logLikelihood <- function(residuals) {
 mean_residual <- mean(residuals)
 sd_residual <- sd(residuals)
 n <- length(residuals)
 
 - n / 2 * log(2 * pi * sd_residual ^ 2) -
  sum((residuals - mean_residual) ^ 2) / (2 * sd_residual ^ 2)
 
}

utils$Likelihood_ratio_index <- function(LogL1, LogL0) {
 #message("Values betweeen 0.2 & 0.4 are indicative of a great fit.")
 1 - (LogL1 / LogL0)
}

utils$Cox_Snell_R2 <- function(LogL1, LogL0, n){
 #message("Maximum value is around 0.75.")
 1 - (exp((LogL0 - LogL1) * (2 / n)))
}
