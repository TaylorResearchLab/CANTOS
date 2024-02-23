compute_silhouette <- function (cluster_df,dist_mat){
  
  cluster_df<-cluster_df[order(cluster_df$SubsetCluster_IDs),]
  dist_mat<-as.data.frame(dist_mat)
  dist_mat$Tumor_Name <- rownames(dist_mat)
  
  for (iter in 1:dim(cluster_df)[1]){
    disease_name <- cluster_df$Tumor_Names[iter]
    cluster_label <- cluster_df$SubsetCluster_IDs[iter]
    cluster_names <- cluster_df$Tumor_Names[which(cluster_df$Tumor_Names %in% )]
    
    subset_distance <- dist_mat %>% dplyr::select()
    
  }
  

  
  return(cluster_df)
}