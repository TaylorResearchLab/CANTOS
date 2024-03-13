suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(jsonlite)
  library(httr)
  library(biomaRt)
  library(ghql)
  library(readxl)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
results_dir <- file.path(analysis_dir,"results")


embedding_v3_large <- read.csv(paste(data_dir,"/embedding_tumor_names_text-embedding-3-large_embeddings.csv",sep=""))
colnames(embedding_v3_large)[1]<-"Tumor_Names"
disease_transform_v3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
disease_transform_v3<-disease_transform_v3[,c(-1)]
rownames(disease_transform_v3)<-embedding_v3_large$Tumor_Name # Needed for AP Clust
rm(embedding_v3_large)


# Set Seed
set.seed(13)
#affinity cluster


#######V3

dist_euclidean_v3<- dist(disease_transform_v3,method = "euclidean")
dist_euclidean_v3<-as.matrix(dist_euclidean_v3)
simmilarity_euclidean_v3<- 1/(1+dist_euclidean_v3)
af_clust_euclidean_v3 <- apcluster(simmilarity_euclidean_v3)#5:11 pm start 7:08 continuing 

affinity_cluster_v3_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_v3_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(af_clust_euclidean_v3@clusters)){
  affinity_cluster_v3_df[iter,1] <- paste(names(unlist(af_clust_euclidean_v3@clusters[iter])),collapse = "@")
  affinity_cluster_v3_df[iter,2] <- iter
}
affinity_cluster_v3_df<- affinity_cluster_v3_df %>% separate_rows(Tumor_Names, sep = '@')
affinity_cluster_v3_df$Cluster_ID<-as.character(affinity_cluster_v3_df$Cluster_ID)
####v3




################ v3 
# Find cluster membership frequencies 
#affinity_cluster_df$Cluster_Total_Members <- NA
cluster_frequency_v3_table <- as.data.frame(table(affinity_cluster_v3_df$Cluster_ID))
colnames(cluster_frequency_v3_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
cluster_frequency_v3_table$Cluster_ID<-as.character(cluster_frequency_v3_table$Cluster_ID)
z_scores_v3<- (cluster_frequency_v3_table$Primary_Cluster_Frequency-mean(cluster_frequency_v3_table$Primary_Cluster_Frequency))/sd(cluster_frequency_v3_table$Primary_Cluster_Frequency)
cluster_frequency_v3_table$z_scores<-z_scores_v3

ind_min_zscore<- which(cluster_frequency_v3_table$z_scores < 2.5)
max_cluster_member <- max(cluster_frequency_v3_table$Primary_Cluster_Frequency[ind_min_zscore])

large_cluster_labels_v3<- cluster_frequency_v3_table$Cluster_ID[which(cluster_frequency_v3_table$Primary_Cluster_Frequency>max_cluster_member)]

converge_list_v3<-list()
disease_transform_v3$Tumor_Names<-rownames(disease_transform_v3)
while(length(large_cluster_labels_v3)>0){
  print(length(large_cluster_labels_v3))
  for(iter in 1:length(large_cluster_labels_v3)){
    Clusters_Names=large_cluster_labels_v3[iter]
    subset_embedding_v3_df <- as.data.frame(affinity_cluster_v3_df$Tumor_Names[affinity_cluster_v3_df$Cluster_ID==Clusters_Names])
    colnames(subset_embedding_v3_df)<-"Tumor_Names"
    rownames(subset_embedding_v3_df)<-subset_embedding_v3_df$Tumor_Names
    subset_embedding_v3_df<- subset_embedding_v3_df %>% dplyr::left_join(disease_transform_v3,by="Tumor_Names")
    rownames(subset_embedding_v3_df)<-subset_embedding_v3_df$Tumor_Names
    subset_embedding_v3_df<-subset_embedding_v3_df[,c(-1)]
    
    result_run_v3_aff<-run_affinity_clustering(Clusters_Names,subset_embedding_v3_df)
    
    flag_converge_v3 <- result_run_v3_aff[[1]]
    subset_affinity_v3_df<-result_run_v3_aff[[2]]
    
    if(flag_converge_v3=="No"){
      for (iter_nested_affinity_cluser_v3 in 1: dim(subset_affinity_v3_df)[1]){
        ind_location_v3 <- which (affinity_cluster_v3_df$Tumor_Names==subset_affinity_v3_df$Tumor_Names[iter_nested_affinity_cluser_v3])
        affinity_cluster_v3_df$Cluster_ID[ind_location_v3]<-subset_affinity_v3_df$SubCluster_ID[iter_nested_affinity_cluser_v3]
      }
    }else if(flag_converge_v3=="Yes"){
      converge_list_v3<-append(Clusters_Names,converge_list_v3)
    }
    
    
  }
  cluster_frequency_v3_table <- as.data.frame(table(affinity_cluster_v3_df$Cluster_ID))
  colnames(cluster_frequency_v3_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
  cluster_frequency_v3_table$Cluster_ID<-as.character(cluster_frequency_v3_table$Cluster_ID)
  large_cluster_labels_v3<- cluster_frequency_v3_table$Cluster_ID[which(cluster_frequency_v3_table$Primary_Cluster_Frequency>max_cluster_member)]
  large_cluster_labels_v3<-setdiff(large_cluster_labels_v3, unlist(converge_list_v3))
}


