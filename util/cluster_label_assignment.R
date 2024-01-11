cluster_label_assignment <- function(affinity_cluster_nested){
  
  
  affinity_cluster_nested$cluster_label <- NA
  subcluster_id_list<- unique(affinity_cluster_nested$SubsetCluster_IDs)
  
  for (iter in 1:length(subcluster_id_list)){
    
    cluster_id_affinity <- subcluster_id_list[iter]
    index_affinity <- which (affinity_cluster_nested$SubsetCluster_IDs==cluster_id_affinity)
    totaling_table <- as.data.frame(table(affinity_cluster_nested$assigned_class[index_affinity]))
    index_max<- which(totaling_table$Freq==max(totaling_table$Freq))
    assigned_cluster_labels <- unique(as.character(totaling_table$Var1[index_max]))
    
    affinity_cluster_nested$cluster_label[index_affinity]<- paste(assigned_cluster_labels,collapse = ";")
  }

  
  
  
  return(affinity_cluster_nested)
}