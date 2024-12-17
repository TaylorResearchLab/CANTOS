compute_embedding_distance <- function(embeddings,distance_name){
  
  distance_matrix<-pdist(embeddings[1:nrow(embeddings),2:4096],metric = distance_name)
  rownames(distance_matrix) <- embeddings$Tumor_Names[1:nrow(embeddings)]
  colnames(distance_matrix) <- embeddings$Tumor_Names[1:nrow(embeddings)]
  return(distance_matrix)
}