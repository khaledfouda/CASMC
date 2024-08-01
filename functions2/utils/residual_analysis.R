utils$residual_analysis <- function(residuals, plot  = TRUE) {
  # input: A vector of residuals
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
