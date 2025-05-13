run_mcnemar_matrix <- function(correct_matrix, adjust_method = "BH") {
  model_names <- colnames(correct_matrix)
  n_models <- ncol(correct_matrix)
  
  # Initialize matrices
  p_matrix <- matrix(NA, nrow = n_models, ncol = n_models,
                     dimnames = list(model_names, model_names))
  n_disagree <- matrix(NA, nrow = n_models, ncol = n_models,
                       dimnames = list(model_names, model_names))
  
  # Store raw p-values in a vector for adjustment
  raw_pvals <- c()
  pair_indices <- list()
  
  for (i in 1:(n_models - 1)) {
    for (j in (i + 1):n_models) {
      model_A <- correct_matrix[, i]
      model_B <- correct_matrix[, j]
      
      tab <- table(factor(model_A, levels = c(1, 0)),
                   factor(model_B, levels = c(1, 0)))
      
      if (all(dim(tab) == c(2, 2))) {
        test <- mcnemar.test(tab, correct = TRUE)
        raw_pvals <- c(raw_pvals, test$p.value)
        pair_indices <- append(pair_indices, list(c(i, j)))
        
        p_matrix[i, j] <- test$p.value
        p_matrix[j, i] <- test$p.value
        
        n_disagree[i, j] <- tab[1, 2] + tab[2, 1]
        n_disagree[j, i] <- n_disagree[i, j]
      }
    }
  }
  
  # Adjust p-values
  adjusted_pvals <- p.adjust(raw_pvals, method = adjust_method)
  
  # Fill adjusted p-value matrix
  p_matrix_adj <- matrix(NA, nrow = n_models, ncol = n_models,
                         dimnames = list(model_names, model_names))
  for (k in seq_along(pair_indices)) {
    i <- pair_indices[[k]][1]
    j <- pair_indices[[k]][2]
    p_matrix_adj[i, j] <- adjusted_pvals[k]
    p_matrix_adj[j, i] <- adjusted_pvals[k]
  }
  diag(p_matrix_adj) <- NA  # Optional: set diagonal to NA
  
  return(list(
    raw_p_values = p_matrix,
    adjusted_p_values = p_matrix_adj,
    disagreements = n_disagree
  ))
}

