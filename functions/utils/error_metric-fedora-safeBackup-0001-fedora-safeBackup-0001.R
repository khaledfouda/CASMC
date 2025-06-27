# generic function:
utils$error_metric = list(
   
#--- error functions:
unexplained_variance = function(predicted, true, adjusted = FALSE, k = NA) {
   # SSE / SST
   if(! adjusted){
      return(sum((true - predicted) ^ 2) / sum((true - mean(true)) ^ 2))
   }else{
      n = length(true)
      stopifnot(is.numeric(k))
      return(((sum((true - predicted) ^ 2) / 
           sum((true - mean(true)) ^ 2)) *
            (n - 1) / (n - k - 1)))
   }
},

mape = function(predicted, true) {
   mean(abs((true - predicted) / true), na.rm = TRUE) * 100
},
mae = function(predicted, true) {
   mean(abs(true - predicted), na.rm = TRUE)
},

rmse_normalized = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE)) / sd(true, na.rm = TRUE)
},

rmse = function(predicted, true) {
   sqrt(mean((true - predicted) ^ 2, na.rm = TRUE))
},

spearman_R2 = function(predicted, true) {
   cor(true, predicted, method = "spearman")
}


)
