run_affinity_clustering <- function(Clusters_Names,subset_embedding_df){
  set.seed(13)
  #affinity_subset <- apcluster(negDistMat(r=2), subset_embedding_df)
  #cat("affinity propogation optimal number of clusters:", length(affinity_subset@clusters), "\n")
 
  dist_euclidean_subset<- dist(subset_embedding_df,method = "euclidean")
  dist_euclidean_subset<-as.matrix(dist_euclidean_subset)
  simmilarity_euclidean_subset<- 1/(1+dist_euclidean_subset)
  affinity_subset <- apcluster(simmilarity_euclidean_subset)
  flag_converge<-"No"
  
  if( length(affinity_subset@clusters) >1 ){
  subset_affinity_df<-as.data.frame(affinity_subset@idx)
  subset_affinity_df<-as.data.frame(matrix(nrow=1,ncol=2))
  colnames(subset_affinity_df)<-c("Tumor_Names","SubCluster_ID")
  
  for (iter_subset in 1: length(affinity_subset@clusters)){
    subset_affinity_df[iter_subset,1] <- paste(names(unlist(affinity_subset@clusters[iter_subset])),collapse = "@")
    subset_affinity_df[iter_subset,2] <- iter_subset
  }
  

  subset_affinity_df<- subset_affinity_df %>% separate_rows(Tumor_Names, sep = '@')
  subset_affinity_df <- subset_affinity_df %>% mutate(SubCluster_ID= paste(Clusters_Names,SubCluster_ID,sep="."))
  }else{
    subset_affinity_df<-as.data.frame(matrix(nrow=dim(subset_embedding_df)[1],ncol=2))
    colnames(subset_affinity_df)<-c("Tumor_Names","SubCluster_ID")
    subset_affinity_df$Tumor_Names<-rownames(subset_embedding_df)
    subset_affinity_df$SubCluster_ID<-Clusters_Names
    flag_converge<-"Yes"
    
  }
  l1<- list(flag_converge,subset_affinity_df)
  return(l1)  

}