residual_analysis <- function(residuals, plot  = TRUE) {
  if (length(residuals) <= 5000)
    shapiro.test(residuals) |> print()
  ks.test(residuals, "pnorm", mean(residuals), sd(residuals)) |> print()
  nortest::ad.test(residuals) |> print()
  
  summary(residuals) |>
    t() |>
    cbind(sd = sd(residuals)) |>
    cbind(count = length(residuals)) |>
    round(5) |>
    print()
  
  if (!plot)
    return()
  
  ggplot(data.frame(residuals), aes(x = residuals)) +
    geom_histogram(bins = 30,
                   fill = "lightblue",
                   color = "black") +
    ggtitle("Histogram of Residuals") +
    xlab("Residuals") +
    ylab("Frequency") -> p1
  
  ggplot(data.frame(residuals), aes(sample = residuals)) +
    geom_qq() +
    geom_qq_line(color = "red") +
    ggtitle("QQ Plot of Residuals") -> p2
  
  ggplot(data.frame(residuals), aes(x = residuals)) +
    geom_density(fill = "lightblue", color = "blue") +
    ggtitle("Density Plot of Residuals") +
    xlab("Residuals") -> p3
  
  ggplot(data.frame(residuals), aes(x = residuals)) +
    geom_boxplot(fill = "lightblue", color = "black") +
    ggtitle("Box Plot of Residuals") +
    xlab("Residuals") -> p4
  
  return(gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2))
}
logLikelihood <- function(residuals) {
  mean_residual <- mean(residuals)
  sd_residual <- sd(residuals)
  n <- length(residuals)
  
  - n / 2 * log(2 * pi * sd_residual ^ 2) -
    sum((residuals - mean_residual) ^ 2) / (2 * sd_residual ^ 2)
  
}

Likelihood_ratio_index <- function(LogL1, LogL0) {
  message("Values betweeen 0.2 & 0.4 are indicative of a great fit.")
  1 - (LogL1 / LogL0)
}

Cox_Snell_R2 <- function(LogL1, LogL0, n){
  message("Maximum value is around 0.75.")
  1 - (exp((LogL0 - LogL1) * (2 / n)))
}
