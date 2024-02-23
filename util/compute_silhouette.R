compute_silhouette <- function (cluster_df,dist_mat){
  
  cluster_df<-cluster_df[order(cluster_df$SubsetCluster_IDs),]
  dist_mat<-as.data.frame(dist_mat)
  dist_mat$Tumor_Name <- rownames(dist_mat)
  cluster_df$silhouette_score<-NA
  
  for (iter in 1:dim(cluster_df)[1]){
    print(iter)
    disease_name <- cluster_df$Tumor_Names[iter]
    cluster_label <- cluster_df$SubsetCluster_IDs[iter]
    cluster_member_names <- cluster_df$Tumor_Names[which(cluster_df$SubsetCluster_IDs==cluster_label)]
      
    subset_dist_mat <- dist_mat %>% dplyr::filter(Tumor_Name %in% cluster_member_names) %>% dplyr::select(any_of(cluster_member_names))
    
    if(dim(subset_dist_mat)[1] > 1){
    
    ind_disease_name<- which(rownames(subset_dist_mat)==disease_name)
    a= subset_dist_mat[ind_disease_name,]
    delete_self_ind <- which(colnames(a)==disease_name)
    a=a[c(-delete_self_ind)]
    a=mean(as.matrix(a))
    
    subset_other_clust <- cluster_df %>% dplyr::filter(SubsetCluster_IDs != cluster_label)
    subset_other_clust<-subset_other_clust[order(subset_other_clust$SubsetCluster_IDs),]
    
    other_clust_dist <- dist_mat %>% dplyr::filter(Tumor_Name %in% disease_name)
    other_clust_dist<- other_clust_dist[,c(-dim(other_clust_dist)[2])]
    other_clust_dist<-t(other_clust_dist)
    colnames(other_clust_dist)<- "dx"
    other_clust_dist<-as.data.frame(other_clust_dist)
    other_clust_dist$dx<-as.double(other_clust_dist$dx)
    other_clust_dist$Tumor_Names <- rownames(other_clust_dist)
    subset_other_clust <- subset_other_clust %>% dplyr::left_join(other_clust_dist,by="Tumor_Names")
    d_clust=aggregate( dx ~ SubsetCluster_IDs,subset_other_clust, mean )
    b= min(d_clust$dx)
    
    silo= (b-a)/max(a,b)
  }else if(dim(subset_dist_mat)[1] == 1){
    silo=0
  }
    cluster_df$silhouette_score[iter] <- silo 
}
  return(cluster_df)
}