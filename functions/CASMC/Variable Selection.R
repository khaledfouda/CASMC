CASMC_var_selection <-
 function(y_train,
          y_valid,
          Y,
          X,
          W_valid,
          error_function = error_metric$rmse,
          # tuning parameters for lambda
          lambda.factor = 1 / 4,
          lambda.init = NULL,
          n.lambda = 20,
          # tuning parameters for J
          rank.init = 2,
          rank.limit = 30,
          rank.step = 2,
          pct = 0.98,
          # laplacian parameters
          lambda.a = 0,
          S.a = NULL,
          lambda.b = 0,
          S.b = NULL,
          # stopping criteria
          early.stopping = 1,
          thresh = 1e-6,
          maxit = 100,
          # trace parameters
          trace = FALSE,
          print.best = TRUE,
          quiet = FALSE,
          # rank constraint parameters
          r_min = 0,
          track = FALSE,
          max_cores = 20,
          # seed
          seed = NULL) {
  k <-  ncol(X)
  correlations <-
   apply(X, 2, function(x)
    mean(abs(cor(
     x, naive_MC(as.matrix(Y))
    )), na.rm = T))
  ordered_indices <- order(correlations, decreasing = TRUE)
  X <- X[, ordered_indices]
  
  remaining_vars <- colnames(X)
  results <-
   data.frame(
    variables = character(),
    var_id = character(),
    validation_error = numeric(),
    rank_beta = numeric(),
    rank_M = numeric(),
    lambda = numeric(),
    J = numeric()
   )
  
  for (i in k:1) {
   X_subset <- X[, 1:i]
   
   tryCatch(
    {
     
   fiti <- CASMC_cv_rank(
    y_train = y_train,
    X = X_subset,
    y_valid = y_valid,
    W_valid = W_valid ,
    trace = trace,
    max_cores = max_cores,
    thresh = thresh,
    lambda.a = lambda.a,
    S.a = S.a,
    lambda.b = lambda.b,
    S.b = S.b,
    lambda.factor = lambda.factor,
    n.lambda = n.lambda,
    rank.limit = rank.limit,
    rank.init  = rank.init,
    rank.step = rank.step,
    pct = pct,
    maxit = maxit,
    r_min = r_min,
    print.best = print.best,
    track  = track,
    seed = seed
   )
   results <- rbind(
    results,
    data.frame(
     rank_beta = fiti$r,
     rank_M = fiti$rank_M,
     lambda = fiti$fit$lambda,
     J = fiti$fit$J,
     validation_error = fiti$error,
     variables = paste(remaining_vars[1:i], collapse =
                        ","),
     var_id = paste(1:i, collapse = ",")
    )
   )
    }, error = function(e) message(paste("Variable Selection (i=",i,"):",e))
   )
  }
  return(results)
 }

remove_collinear_cols <- function(matrix, thresh=0.7){
  
  cor_m <- cor(matrix)
  to_remove <- caret::findCorrelation(cor_m, thresh, verbose=TRUE)  
  matrix <- matrix[, - to_remove, drop=FALSE]
  return(matrix)
}













