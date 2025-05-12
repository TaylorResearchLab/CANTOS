run_mcnemar_matrix <- function(correct_matrix) {
  model_names <- colnames(correct_matrix)
  n_models <- ncol(correct_matrix)
  
  # Initialize matrices to store p-values and disagreement counts
  p_matrix <- matrix(NA, nrow = n_models, ncol = n_models,
                     dimnames = list(model_names, model_names))
  n_disagree <- matrix(NA, nrow = n_models, ncol = n_models,
                       dimnames = list(model_names, model_names))
  
  for (i in 1:(n_models - 1)) {
    for (j in (i + 1):n_models) {
      model_A <- correct_matrix[, i]
      model_B <- correct_matrix[, j]
      
      # Build 2x2 contingency table
      tab <- table(factor(model_A, levels = c(1, 0)),
                   factor(model_B, levels = c(1, 0)))
      
      # Handle rare case of missing categories
      if (all(dim(tab) == c(2, 2))) {
        test <- mcnemar.test(tab, correct = TRUE)
        p_matrix[i, j] <- test$p.value
        p_matrix[j, i] <- test$p.value
        
        n_disagree[i, j] <- tab[1, 2] + tab[2, 1]  # off-diagonal disagreements
        n_disagree[j, i] <- n_disagree[i, j]
      }
    }
  }
  
  return(list(p_values = p_matrix, disagreements = n_disagree))
}