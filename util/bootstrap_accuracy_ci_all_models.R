bootstrap_accuracy_ci_all_models <- function(correct_matrix, n_iter = 10000, conf_level = 0.95) {
  n_models <- ncol(correct_matrix)
  model_names <- colnames(correct_matrix)
  n_samples <- nrow(correct_matrix)
  
  # Initialize results
  results <- data.frame(
    Model = model_names,
    Accuracy = NA_real_,
    CI_Lower = NA_real_,
    CI_Upper = NA_real_,
    stringsAsFactors = FALSE
  )
  
  alpha <- 1 - conf_level
  
  for (i in 1:n_models) {
    accs <- numeric(n_iter)
    
    # Bootstrap loop
    for (j in 1:n_iter) {
      idx <- sample(1:n_samples, replace = TRUE)
      accs[j] <- mean(correct_matrix[idx, i])
    }
    
    # Store mean and CI
    results$Accuracy[i]  <- mean(accs)
    results$CI_Lower[i]  <- quantile(accs, probs = alpha / 2)
    results$CI_Upper[i]  <- quantile(accs, probs = 1 - alpha / 2)
  }
  
  # Optional: round for clean reporting
  results$Accuracy <- round(results$Accuracy, 3)
  results$CI_Lower <- round(results$CI_Lower, 3)
  results$CI_Upper <- round(results$CI_Upper, 3)
  
  return(results)
}
