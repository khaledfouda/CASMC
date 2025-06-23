utils$MC_train_test_split <-
   function(obs_mask,
            testp = 0.2,
            seed = NULL) {
      #' returns a new mask similar to mask with a new train-test sets.
      #' testp: proportion of test out of the nonzero cells in mask
      #' new_mask:
      #------------------------------------------------------------------
      #  1 -> train    (obs_mask=1) |
      #  0 -> test     (obs_mask=1) | 
      #  1 -> missing  (obs_mask=0) 
      #-----------------------------------------------------------------
      if(!is.null(seed)) set.seed(seed)
      n_rows <- dim(obs_mask)[1]
      n_cols <- dim(obs_mask)[2]
      # Create a data frame of all matrix indices
      # we only consider non-missing data (ie, with mask_ij=1)
      indices <- expand.grid(row = 1:n_rows, col = 1:n_cols)[obs_mask == 1, ]
      # Shuffle indices (both rows and columns are shuffled. later, we will reshuffle the columns)
      indices <-  indices[sample(1:nrow(indices)),]
      row.names(indices) <- NULL
      
      test.idx = sample(1:nrow(indices),
                        size = nrow(indices) * testp,
                        replace = FALSE)
      test.indices = indices[test.idx, ]
      
      new_mask = matrix(1, nrow = n_rows, ncol = n_cols)
      new_mask[obs_mask == 0] <- 1
      new_mask[as.matrix(test.indices[, c("row", "col")])] <- 0
      
      return(new_mask)
   }

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
         #  1 -> train    (obs_mask=1) |
         #  0 -> valid    (obs_mask=1) | 
         #  1 -> missing  (obs_mask=0) (missing)
         valid_mask <- matrix(1, nrow = n_rows, ncol = n_cols)
         test_indices <- indices[indices$fold == i,]
         valid_mask[obs_mask == 0] <- 1
         valid_mask[as.matrix(test_indices[, c("row", "col")])] <- 0
         # Store the mask
         folds[[i]] <- valid_mask
      }
      return(folds)
   }

