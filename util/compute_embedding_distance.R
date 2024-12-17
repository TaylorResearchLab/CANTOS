compute_embedding_distance <- function(embeddings,distance_name){
  
  distance_matrix<-pdist(embeddings,metric = distance_name)
  return(distance_matrix)
}