nested_clust_edit_dist <- function(n,affinity_cluster_df,dist_mat){

  affinity_cluster_df$SubsetCluster_IDs<-affinity_cluster_df$Cluster_ID
  dist_mat$Tumor_Name <- rownames(dist_mat)
  
  for(iter_nest in 1:n){
    affinity_cluster_df$SubsetCluster_IDs <- as.character(affinity_cluster_df$SubsetCluster_IDs)
    affinity_subcluster_labels <- unique(affinity_cluster_df$SubsetCluster_IDs)
    
    for(iter in 1:length(affinity_subcluster_labels)){
      nested_dist_df <- as.data.frame(affinity_cluster_df$Tumor_Names[affinity_cluster_df$SubsetCluster_IDs==affinity_subcluster_labels[iter]])
      
      colnames(nested_dist_df)<-"Tumor_Name"
      rownames(nested_dist_df)<-nested_dist_df$Tumor_Name
      nested_dist_df<- nested_dist_df %>% dplyr::left_join(dist_mat,by="Tumor_Name")
      rownames(nested_dist_df)<-nested_dist_df$Tumor_Name
      nested_dist_df <- nested_dist_df %>% dplyr::select(any_of(rownames(nested_dist_df)))

      if(dim(nested_dist_df)[1]>2){
        
        affinity_subset <- apcluster(nested_dist_df)
        cat("affinity propogation optimal number of clusters:", length(affinity_subset@clusters), "\n")
        
        nested_subset_affinity_df<-as.data.frame(affinity_subset@idx)
        nested_subset_affinity_df<-as.data.frame(matrix(nrow=1,ncol=2))
        colnames(nested_subset_affinity_df)<-c("Tumor_Names","SubCluster_ID")
        for (iter_subset in 1: length(affinity_subset@clusters)){
          nested_subset_affinity_df[iter_subset,1] <- paste(names(unlist(affinity_subset@clusters[iter_subset])),collapse = "@")
          nested_subset_affinity_df[iter_subset,2] <- iter_subset
        }
        nested_subset_affinity_df<- nested_subset_affinity_df %>% separate_rows(Tumor_Names, sep = '@')
        nested_subset_affinity_df <- nested_subset_affinity_df %>% mutate(SubCluster_ID= paste(affinity_subcluster_labels[iter],SubCluster_ID,sep="."))
        
        for (iter_nested_affinity_cluser in 1: dim(nested_subset_affinity_df)[1]){
          ind_location <- which (affinity_cluster_nested$Tumor_Names==nested_subset_affinity_df$Tumor_Names[iter_nested_affinity_cluser])
          affinity_cluster_nested$SubsetCluster_IDs[ind_location]<-nested_subset_affinity_df$SubCluster_ID[iter_nested_affinity_cluser]
        }
      }
      
    }
    
  }
  
  return(affinity_cluster_df)
}