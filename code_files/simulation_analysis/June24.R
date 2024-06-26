dat <-
 generate_simulation_rows(
  600,
  700,
  r = 10,
  k = 10, 
  missing_prob = 0.9,
  coll = T,
  prepare_for_fitting = TRUE,
  half_discrete = FALSE,
  informative_cov_prop = 1,mar_sparse = T,
  mv_beta = T,
  seed = 2023
 )