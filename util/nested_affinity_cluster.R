nested_affinity_clusted <- function(n,affinity_cluster_nested){
  
  for (iter_nest in 1:n){
    
    affinity_cluster_nested$SubsetCluster_IDs <- as.character(affinity_cluster_nested$SubsetCluster_IDs)
    
    affinity_subcluster_labels <- unique(affinity_cluster_nested$SubsetCluster_IDs)
    
    for(iter in 1:length(affinity_subcluster_labels)){
      
      nested_embedding_df <- as.data.frame(affinity_cluster_nested$Tumor_Names[affinity_cluster_nested$SubsetCluster_IDs==affinity_subcluster_labels[iter]])
      
      colnames(nested_embedding_df)<-"Tumor_Name"
      rownames(nested_embedding_df)<-nested_embedding_df$Tumor_Name
      nested_embedding_df<- nested_embedding_df %>% dplyr::left_join(disease_transform,by="Tumor_Name")
      rownames(nested_embedding_df)<-nested_embedding_df$Tumor_Name
      nested_embedding_df<-nested_embedding_df[,c(-1)]
      
      if(dim(nested_embedding_df)[1]>2){
        
        affinity_subset <- apcluster(negDistMat(r=2), nested_embedding_df)
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
  return(affinity_cluster_nested)
  
}