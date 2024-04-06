cluster_label_assignment <- function(affinity_cluster_nested){
  
  
  affinity_cluster_nested$who_cluster_label<- NA
  affinity_cluster_nested$ncit_cluster_label<- NA
  
  
  subcluster_id_list<- unique(affinity_cluster_nested$Cluster_ID)
  
  for (iter in 1:length(subcluster_id_list)){
    
    cluster_id_affinity <- subcluster_id_list[iter]
    index_affinity <- which (affinity_cluster_nested$Cluster_ID==cluster_id_affinity)
    #totaling_table <- as.data.frame(table(affinity_cluster_nested$assigned_class[index_affinity]))
    
    totaling_table_who <- as.data.frame(table(affinity_cluster_nested$WHO_Matches[index_affinity]))
    totaling_table_ncit <- as.data.frame(table(affinity_cluster_nested$NCIT_Matches[index_affinity]))
    
    
    index_max_who<- which(totaling_table_who$Freq==max(totaling_table_who$Freq))
    index_max_ncit<- which(totaling_table_ncit$Freq==max(totaling_table_ncit$Freq))
    
    if(length(index_max_who)==1){
      affinity_cluster_nested$who_cluster_label[index_affinity]<- as.character(totaling_table_who$Var1[index_max_who])
    }else{
      for(iter_who in 1:length(index_max_who)){
        if(affinity_cluster_nested$NCIT_Matches)
      }
    }
    
    if(length(index_max_ncit)==1){
      affinity_cluster_nested$ncit_cluster_label[index_affinity]<- as.character(totaling_table_who$Var1[index_max_ncit])
    }else{
      
    }
    
    who_cluster_labels <- unique(as.character(totaling_table_who$Var1[index_max_who]))
    ncit_cluster_labels <- unique(as.character(totaling_table_ncit$Var1[index_max_ncit]))
    
    affinity_cluster_nested$who_cluster_label[index_affinity]<- paste(who_cluster_labels,collapse = ";")
    affinity_cluster_nested$ncit_cluster_label[index_affinity]<- paste(ncit_cluster_labels,collapse = ";")
    
  }
  
  
  
  
  return(affinity_cluster_nested)
}