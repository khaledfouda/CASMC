# Computing the test error as given by Mao in page 205
unexplained_variance <- function(predicted, actual){
  sum( (actual-predicted)^2 )/ sum((actual-mean(actual))^2)
}


MAPE_error <- function(predicted, actual) {
    mean(abs((actual - predicted) / actual), na.rm = TRUE) * 100
}


adjusted_unexplained_variance <- function(predicted, actual, p=1, n=length(actual)) {
   
      ((sum( (actual-predicted)^2 )/ sum((actual-mean(actual))^2)) * (n - 1) / (n - p - 1))
}

normalized_RMSE <- function(predicted, actual) {
   sqrt(mean((actual - predicted)^2, na.rm = TRUE)) / sd(actual, na.rm = TRUE)
}

RMSE_error <- function(predicted, actual) {
   sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

mao_error <- function(predicted, actual){
   # the test error defined by Mao in the simulation section
   sqrt( sum((predicted-actual)^2) / sum(actual^2))
}


# watch out! functions send predicted, then actual
# Last line is run last
#test_error <- adjusted_unexplained_variance
test_error <<- RMSE_error
#test_error <- normalized_RMSE
#test_error <- MAPE

