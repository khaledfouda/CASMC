
utils$MC_Kfold_split <-
 function(n_rows,
          n_cols,
          n_folds,
          obs_mask,
          seed = NULL) {
  
  if(! is.null(seed)) set.seed(seed)
  # Create a data frame of all matrix indices
  indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)
  # we only consider non-missing data (ie, with obs_mask=1)
  indices <- indices[obs_mask == 1, ]
  # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
  indices <- indices[sample(1:nrow(indices)),]
  
  # Assign each observed index to one of k groups, 
  # ensuring that the number of validation cells in each row are equal (or close)
  indices <- indices %>%
   group_by(row) %>%
   # the following is to shuffle within each row
   do(sample_n(., size = nrow(.))) %>%
   mutate(fold = rep(1:n_folds, length.out = n())) %>%
   ungroup()
  
  # Assign each index to one of k groups
  # Create a list to hold each fold
  folds <- vector("list", n_folds)
  for (i in 1:n_folds) {
   # Create a mask for the test cells in this fold
   #  1 -> train (obs_mask=1) |
   #  0 -> valid (obs_mask=1) | 
   #  0 -> missing  (obs_mask=0) (missing)
   valid_mask <- matrix(1, nrow = n_rows, ncol = n_cols)
   test_indices <- indices[indices$fold == i,]
   valid_mask[obs_mask == 0] <- 0
   valid_mask[as.matrix(test_indices[, c("row", "col")])] <- 0
   # Store the mask
   folds[[i]] <- valid_mask
  }
  return(folds)
 }

