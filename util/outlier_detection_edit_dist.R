outlier_detection_edit_dist <- function(nested_affinity_cluster,dissimilarity_matrix){
  
  cluster_label <- unique(nested_affinity_cluster$Cluster_ID)
  dissimilarity_matrix<- as.data.frame(dissimilarity_matrix)
  nested_affinity_cluster$isolation_outlier_score<-NA
  for(iter in 1:length(cluster_label)){
    current_cluster_label <- cluster_label[iter]
    subset_df <- nested_affinity_cluster %>% filter(Cluster_ID==current_cluster_label) %>% dplyr::select(Tumor_Names,Cluster_ID)
    if(dim(subset_df)[1]>2){
      subset_distance <- dissimilarity_matrix %>% dplyr::select(one_of(subset_df$Tumor_Names))
      subset_distance$Tumor_Names<-rownames(subset_distance)
      subset_distance <- subset_distance %>% filter(Tumor_Names %in% colnames(subset_distance))
      subset_df<-subset_df %>% dplyr::left_join(subset_distance,by="Tumor_Names")
      
      #model <- isolation.forest(subset_df[1:nrow(subset_df),3:ncol(subset_df)], ndim=3, ntrees=ceiling(sqrt(ncol(subset_df)-2)), nthreads=1) # ntrees 50 initially
      
      model <- isolation.forest(subset_df[1:nrow(subset_df),3:ncol(subset_df)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
      
      scores <- predict(model, subset_df[1:nrow(subset_df),3:ncol(subset_df)], type="score")
      ind_clust <- which(nested_affinity_cluster$Cluster_ID==current_cluster_label)
      nested_affinity_cluster$isolation_outlier_score[ind_clust]<-scores
    }else{
      ind_clust <- which(nested_affinity_cluster$Cluster_ID==current_cluster_label)
      nested_affinity_cluster$isolation_outlier_score[ind_clust]<-0
    }
    
    
  }
  nested_affinity_cluster<- nested_affinity_cluster %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))
  
  ###### LOF Computation
  
  nested_affinity_cluster$LOF_Scores<-NA
  lof_scores_minpts_list<-list()
  
  for(iter in 1:length(cluster_label)){
    current_cluster_label <- cluster_label[iter]
    ind_clust <- which(nested_affinity_cluster$Cluster_ID==current_cluster_label)
    lof_scores_minpts_list<-list()
    
    subset_df <- nested_affinity_cluster %>% filter(Cluster_ID==current_cluster_label) %>% dplyr::select(Tumor_Names,Cluster_ID)
    
    subset_distance <- dissimilarity_matrix %>% dplyr::select(one_of(subset_df$Tumor_Names))
    subset_distance$Tumor_Names<-rownames(subset_distance)
    subset_distance <- subset_distance %>% filter(Tumor_Names %in% colnames(subset_distance))
    subset_df<-subset_df %>% dplyr::left_join(subset_distance,by="Tumor_Names")
    
    
    if(dim(subset_df)[1]>2){ # Need at least 2 data points to run isolation forest
      
      min_pts<- 2:(dim(subset_df)[1]-1)
      for(iter_pts in min_pts){
        lof_scores_minpts <- lof(subset_df[,3:ncol(subset_df)],iter_pts)
        lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
      }
      lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
      lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
      nested_affinity_cluster$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
      
    }else{
      nested_affinity_cluster$LOF_Scores[ind_clust]<-0
    }
    
    
  }
  nested_affinity_cluster<- nested_affinity_cluster %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))
  
  return(nested_affinity_cluster)
  
}
  