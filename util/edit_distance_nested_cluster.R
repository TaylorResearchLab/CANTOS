edit_distance_nested_cluster <- function(affinity_cluster_df,dist_mat){
  
  #affinity_cluster_df$SubsetCluster_IDs<-affinity_cluster_df$Cluster_ID
  dist_mat<-as.data.frame(dist_mat)
  dist_mat$Tumor_Names <- rownames(dist_mat)
  
  
  cluster_frequency_table<- as.data.frame(table(affinity_cluster_df$Cluster_ID))
  colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
  cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
  
  max_cluster_member<- median(cluster_frequency_table$Primary_Cluster_Frequency)
  
  large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]
 
  while(length(large_cluster_labels)>0 & flag !=0 ){
  print(length(large_cluster_labels))
  for (iter in 1:length(large_cluster_labels)){
    current_cluster_label <- large_cluster_labels[iter]
    nested_cluster_df<- affinity_cluster_df %>% dplyr::filter(Cluster_ID==current_cluster_label)%>%dplyr::select(Tumor_Names)
    nested_cluster_df <- nested_cluster_df %>% left_join(dist_mat,by="Tumor_Names")
    
    rownames(nested_cluster_df)<-nested_cluster_df$Tumor_Names
    nested_cluster_df <- nested_cluster_df %>% dplyr::select(any_of(rownames(nested_cluster_df)))
    rownames(nested_cluster_df)<-colnames(nested_cluster_df)
    nested_cluster_df<-as.matrix(nested_cluster_df)
    
    if(dim(nested_cluster_df)[1]>2){
      
      affinity_subset <- apcluster(nested_cluster_df)
      nested_subset_affinity_df<-as.data.frame(affinity_subset@idx)
      nested_subset_affinity_df<-as.data.frame(matrix(nrow=1,ncol=2))
      colnames(nested_subset_affinity_df)<-c("Tumor_Names","Cluster_ID")
      for (iter_subset in 1: length(affinity_subset@clusters)){
        nested_subset_affinity_df[iter_subset,1] <- paste(names(unlist(affinity_subset@clusters[iter_subset])),collapse = "@")
        nested_subset_affinity_df[iter_subset,2] <- iter_subset
      }
      
      nested_subset_affinity_df<- nested_subset_affinity_df %>% separate_rows(Tumor_Names, sep = '@')
      nested_subset_affinity_df <- nested_subset_affinity_df %>% mutate(Cluster_ID= paste(current_cluster_label,Cluster_ID,sep="."))
      
      for (iter_nested_affinity_cluser in 1: dim(nested_subset_affinity_df)[1]){
        ind_location <- which (affinity_cluster_df$Tumor_Names==nested_subset_affinity_df$Tumor_Names[iter_nested_affinity_cluser])
        affinity_cluster_df$Cluster_ID[ind_location]<-nested_subset_affinity_df$Cluster_ID[iter_nested_affinity_cluser]
      }
    
  }
  
  
  }
    cluster_frequency_table<- as.data.frame(table(affinity_cluster_df$Cluster_ID))
    colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
    cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
    large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]

}
  
  return(affinity_cluster_df)
  
}