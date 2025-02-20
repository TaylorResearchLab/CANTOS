compute_mi_pairwise<-function(data_mat){
  
  num_methods<-ncol(data_mat)
  MI_mat <-matrix(0, nrow = num_methods, ncol = num_methods)
  for (i in 1:num_methods) {
    for (j in i:num_methods) {
      mi_score <- mutinformation(data_mat[[i]], data_mat[[j]])
      MI_mat[i, j] <- mi_score
      MI_mat[j, i] <- mi_score  # Fill both halves since it's symmetric
    }
  }
  colnames(MI_mat) <- colnames(data_mat)
  rownames(MI_mat) <- colnames(data_mat)
  return(MI_mat)
}