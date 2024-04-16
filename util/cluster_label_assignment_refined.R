cluster_label_assignment_refined <- function(affinity_cluster_nested){
  
  
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
         high_freq_labels<-as.character(totaling_table_who$Var1[index_max_who])
         other_freq_label <- affinity_cluster_nested %>% filter(Cluster_ID==cluster_id_affinity)%>%dplyr::select(WHO_Matches)%>%distinct()
         other_freq_label<-other_freq_label$WHO_Matches
         other_freq_label<-setdiff(other_freq_label,high_freq_labels)
         
         for(iter_high_freq in 1:length(high_freq_labels)){
           index_locator_high<-which(affinity_cluster_nested$WHO_Matches==high_freq_labels[iter_high_freq])
           index_locator_high <- intersect(index_locator_high,index_affinity)
           affinity_cluster_nested$who_cluster_label[index_locator_high]<-high_freq_labels[iter_high_freq]
         }
         if(length(other_freq_label)>=1){
           for(iter_other_freq in 1:length(other_freq_label)){
             index_locator_other<-which(affinity_cluster_nested$WHO_Matches==other_freq_label[iter_other_freq])
             index_locator_other <- intersect(index_locator_other,index_affinity)
             affinity_cluster_nested$who_cluster_label[index_locator_other]<-other_freq_label[iter_other_freq]
             
           }
         }
         
      }
    
    
    if(length(index_max_ncit)==1){
      affinity_cluster_nested$ncit_cluster_label[index_affinity]<- as.character(totaling_table_ncit$Var1[index_max_ncit])
    }else{
      high_freq_labels_ncit<-as.character(totaling_table_ncit$Var1[index_max_ncit])
      other_freq_label_ncit <- affinity_cluster_nested %>% filter(Cluster_ID==cluster_id_affinity)%>%dplyr::select(NCIT_Matches)%>%distinct()
      other_freq_label_ncit<-other_freq_label_ncit$NCIT_Matches
      other_freq_label_ncit<-setdiff(other_freq_label_ncit,high_freq_labels_ncit)
      
      for(iter_high_freq in 1:length(high_freq_labels_ncit)){
        index_locator_high<-which(affinity_cluster_nested$NCIT_Matches==high_freq_labels_ncit[iter_high_freq])
        index_locator_high <- intersect(index_locator_high,index_affinity)
        affinity_cluster_nested$ncit_cluster_label[index_locator_high]<-high_freq_labels_ncit[iter_high_freq]
      }
      if(length(other_freq_label_ncit)>=1){
        for(iter_other_freq in 1:length(other_freq_label_ncit)){
          index_locator_other<-which(affinity_cluster_nested$NCIT_Matches==other_freq_label_ncit[iter_other_freq])
          index_locator_other <- intersect(index_locator_other,index_affinity)
          affinity_cluster_nested$ncit_cluster_label[index_locator_other]<-other_freq_label_ncit[iter_other_freq]
          
        }
      }
      
    }
    
 
    
  }
  
  
  
  
  return(affinity_cluster_nested)
}