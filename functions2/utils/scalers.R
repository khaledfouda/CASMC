


utils$scalers <-
 list(
  median = function(X, skip_bin = TRUE) {
   fun <- function(x) {
    (x - median(x)) / IQR(x)
   }
   utils$scalers_post_process(X, fun, skip_bin)
  },
  
  minmax = function(X, skip_bin = TRUE) {
   fun <- function(x) {
    (x - min(x)) / (max(x) - min(x))
   }
   utils$scalers_post_process(X, fun, skip_bin)
  },
  
  standard = function(X, skip_bin = TRUE) {
   fun <- function(x) {
    (x - mean(x)) / sd(x)
   }
   utils$scalers_post_process(X, fun, skip_bin)
  },
  
  positive_log1p = function(X, skip_bin = TRUE) {
   fun <- function(x) {
    log1p(x + ifelse(min(x) < 0, abs(min(x)), 0))
   }
   utils$scalers_post_process(X, fun, skip_bin)
  }
 )



utils$scalers_post_process <- function(X, fun, skip_bin) {
 is.mat = FALSE
 if (!is.null(ncol(X)) && ncol(X) > 1)
  is.mat = TRUE
 
 if (is.mat == FALSE) {
  if (length(unique(X)) <= 2 && skip_bin) {
   return(X)
  } else
   return(fun(X))
 } else{
  fun2 = function(x, f, skip_bin) {
   if (length(unique(x)) <= 2 && skip_bin) {
    return(x)
   } else
    return(fun(x))
  }
  return(apply(X, 2, fun2, f = fun, skip_bin = skip_bin))
  
 }
 
}