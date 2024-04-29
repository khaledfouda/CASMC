# generic function:
error_metric = list(
   
#--- error functions:
unexplained_variance = function(predicted, true) {
   sum((true - predicted) ^ 2) / sum((true - mean(true)) ^ 2)
},

adjusted_unexplained_variance = function(predicted, true, p = 1, n = length(true)) {
   ((sum((true - predicted) ^ 2) / 
        sum((true - mean(true)) ^ 2)) *
       (n - 1) / (n - p - 1))
},

mape = function(predicted, true) {
   mean(abs((true - predicted) / true), na.rm = TRUE) * 100
},

rmse_normalized = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE)) / sd(true, na.rm = TRUE)
},

rmse = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE))
},

rmse_mao = function(predicted, true) {
   # a form of normalized RMSE that was used in Mao's simulation
   sqrt(sum((predicted - true) ^ 2) / sum(true ^ 2))
},

spearman = function(predicted, true) {
   cor(true, predicted, method = "spearman")
}


)
