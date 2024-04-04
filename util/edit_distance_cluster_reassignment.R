edit_distance_cluster_reassignment <- function (nested_affinity_cluster){
  nested_affinity_cluster$Cluster_ID<-as.character(nested_affinity_cluster$Cluster_ID)
  
  nested_affinity_cluster_reassigned_df<-nested_affinity_cluster %>% dplyr::select(Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)
  
  outlier_cluster <- nested_affinity_cluster%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Cluster_ID)%>%distinct()
  
  
  
  for(iter in 1:dim(outlier_cluster)[1]){
    current_cluster_id <- outlier_cluster$Cluster_ID[iter]
    
    cluster_subset <- nested_affinity_cluster%>%filter(Cluster_ID==current_cluster_id)
    
    ind_outliers<- which(cluster_subset$Isolation_Outlier=="Yes" | cluster_subset$LOF_Outlier=="Yes")
    
    new_suffix_cluster_id <-1:length(ind_outliers)
    new_cluster_id<- paste(current_cluster_id,new_suffix_cluster_id,sep=";")
    cluster_subset$Cluster_ID[ind_outliers]<-new_cluster_id
    
    cluster_subset<-cluster_subset%>%dplyr::select(Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)
    
    table_frequency_cluster_id<- as.data.frame(table(cluster_subset$Cluster_ID))
    table_frequency_cluster_id$Var1<-as.character(table_frequency_cluster_id$Var1)
    high_frequency_cluster_id<- table_frequency_cluster_id$Var1[which(table_frequency_cluster_id$Freq>1)]
    
    #low_frequency_cluster_id<-table_frequency_cluster_id$Var1[which(table_frequency_cluster_id$Freq==1)]
    #low_frequency_assigned_class <- cluster_subset %>% filter(Cluster_ID %in% low_frequency_cluster_id) %>% dplyr::select(assigned_class)
    
    for(iter_subset in 1:dim(cluster_subset)[1]){
      if(!(cluster_subset$Cluster_ID[iter_subset] %in% high_frequency_cluster_id)){
        cluster_subset$who_cluster_label[iter_subset]<-cluster_subset$WHO_Matches[iter_subset]
        cluster_subset$ncit_cluster_label[iter_subset]<-cluster_subset$NCIT_Matches[iter_subset]
        
      }else{
        indx_ammend<- which(cluster_subset$Cluster_ID==cluster_subset$Cluster_ID[iter_subset])
        all_who_class<- unique(cluster_subset$WHO_Matches[indx_ammend])
        all_who_class<-paste(all_who_class,collapse=";")
        cluster_subset$who_cluster_label[iter_subset]<-all_who_class
        
        all_ncit_class<- unique(cluster_subset$NCIT_Matches[indx_ammend])
        all_ncit_class<-paste(all_ncit_class,collapse=";")
        cluster_subset$ncit_cluster_label[iter_subset]<-all_ncit_class
      }
    }
    
    
    
    nested_affinity_cluster_reassigned_df<-nested_affinity_cluster_reassigned_df %>% rows_update(cluster_subset,by="Tumor_Names")
  }
  
  cluster_labels<-unique(nested_affinity_cluster_reassigned_df$Cluster_ID)
  cluster_labels<- as.data.frame(cbind(cluster_labels,c(1:length(cluster_labels))))
  colnames(cluster_labels)<-c("Cluster_ID","updated_ID")
  nested_affinity_cluster_reassigned_df<-nested_affinity_cluster_reassigned_df %>% dplyr::left_join(cluster_labels,by="Cluster_ID")
  nested_affinity_cluster_reassigned_df<-nested_affinity_cluster_reassigned_df[,c(1,9,3:8,2)]
  
  nested_affinity_cluster_reassigned_df$updated_ID<-as.numeric(nested_affinity_cluster_reassigned_df$updated_ID)
  
  return(nested_affinity_cluster_reassigned_df)
  
  
}