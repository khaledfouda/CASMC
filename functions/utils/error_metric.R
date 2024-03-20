# Computing the test error as given by Mao in page 205
unexplained_variance <- function(predicted, true) {
   sum((true - predicted) ^ 2) / sum((true - mean(true)) ^ 2)
}


MAPE_error <- function(predicted, true) {
   mean(abs((true - predicted) / true), na.rm = TRUE) * 100
}


adjusted_unexplained_variance <-
   function(predicted,
            true,
            p = 1,
            n = length(true)) {
      ((sum((true - predicted) ^ 2) / sum((true - mean(
         true
      )) ^ 2)) * (n - 1) / (n - p - 1))
   }

normalized_RMSE <- function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE)) / sd(true, na.rm = TRUE)
}

RMSE_error <- function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE))
}

mao_error <- function(predicted, true) {
   # the test error defined by Mao in the simulation section
   sqrt(sum((predicted - true) ^ 2) / sum(true ^ 2))
}


# watch out! functions send predicted, then true
# Last line is run last
#test_error <- adjusted_unexplained_variance
test_error <<- RMSE_error
#test_error <- normalized_RMSE
#test_error <- MAPE
