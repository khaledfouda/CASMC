CASMC_Lasso_hparams <-
  list(
    M = list(
      # tuning parameters for lambda
      lambda.factor = 1 / 4,
      lambda.init = NULL,
      n.lambda = 20,
      # tuning parameters for J
      rank.init = 2,
      rank.limit = 30,
      rank.step = 2,
      pct = 0.98,
      early.stopping = 1
    ),
    beta = list(
      # L1 parameters
      learning.rate = "default",
      lambda.max = NULL,
      prox.iter.max = 20,
      n.lambda = 20
    ),
    laplacian = list(
      # laplacian parameters
      lambda.a = 0,
      S.a = NULL,
      lambda.b = 0,
      S.b = NULL
    )
  )