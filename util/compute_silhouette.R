compute_silhouette <- function (cluster_df,dist_mat){
  
  cluster_df<-cluster_df[order(cluster_df$SubsetCluster_IDs),]
  dist_mat<-as.data.frame(dist_mat)
  dist_mat$Tumor_Name <- rownames(dist_mat)
  
  for (iter in 1:dim(cluster_df)[1]){
    disease_name <- cluster_df$Tumor_Names[iter]
    cluster_label <- cluster_df$SubsetCluster_IDs[iter]
    cluster_member_names <- cluster_df$Tumor_Names[which(cluster_df$SubsetCluster_IDs==cluster_label)]
      
    subset_dist_mat <- dist_mat %>% dplyr::filter(Tumor_Name %in% cluster_member_names) %>% dplyr::select(any_of(cluster_member_names))
    
    ind_disease_name<- which(rownames(subset_dist_mat)==disease_name)
    a= subset_dist_mat[ind_disease_name,]
    delete_self_ind <- which(colnames(a)==disease_name)
    a=a[c(-delete_self_ind)]
    a=mean(as.matrix(a))
    
    subset_other_clust <- cluster_df %>% dplyr::filter(SubsetCluster_IDs != cluster_label)
    subset_other_clust<-subset_other_clust[order(subset_other_clust$SubsetCluster_IDs),]
    
    

  }
  

  
  return(cluster_df)
}