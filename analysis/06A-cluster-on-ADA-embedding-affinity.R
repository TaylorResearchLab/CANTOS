# Computes affinity cluster of ADA2 data. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.
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
  library(apcluster)
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
source(paste(util_dir,"/run_affinity_clustering.R",sep=""))

########################################*************************************###########
# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2.csv",sep="") )
colnames(disease_transform)[1]<-"Tumor_Name"
rownames(disease_transform)<-disease_transform$Tumor_Name # Needed for AP Clust

# Set Seed
set.seed(13)
#affinity cluster

dist_euclidean<- dist(disease_transform,method = "euclidean")
dist_euclidean<-as.matrix(dist_euclidean)
simmilarity_euclidean<- 1/(1+dist_euclidean)
af_clust_euclidean <- apcluster(simmilarity_euclidean) # 2:01 pm
cat("affinity propogation optimal number of clusters:", length(af_clust_euclidean@clusters), "\n")#3550 clusters

#d.apclus2 <- apcluster(negDistMat(r=2), disease_transform) # 1 hr 28 mins 11:28 pm - 12:08 pm
#cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n") #1113 clusters 
affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(af_clust_euclidean@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(af_clust_euclidean@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
affinity_cluster_df$Cluster_ID<-as.character(affinity_cluster_df$Cluster_ID)

# Find cluster membership frequencies 
#affinity_cluster_df$Cluster_Total_Members <- NA
cluster_frequency_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
z_scores<- (cluster_frequency_table$Primary_Cluster_Frequency-mean(cluster_frequency_table$Primary_Cluster_Frequency))/sd(cluster_frequency_table$Primary_Cluster_Frequency)
cluster_frequency_table$z_scores<-z_scores

ind_min_zscore<- which(cluster_frequency_table$z_scores < 2.5)
max_cluster_member <- max(cluster_frequency_table$Primary_Cluster_Frequency[ind_min_zscore])

#median_cluster_frequency <- median(cluster_frequency_table$Primary_Cluster_Frequency)

#large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>median_cluster_frequency)]
large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]

converge_list<-list()

while(length(large_cluster_labels)>0){
  print(length(large_cluster_labels))
  for(iter in 1:length(large_cluster_labels)){
    Clusters_Names=large_cluster_labels[iter]
    subset_embedding_df <- as.data.frame(affinity_cluster_df$Tumor_Names[affinity_cluster_df$Cluster_ID==Clusters_Names])
    colnames(subset_embedding_df)<-"Tumor_Name"
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<- subset_embedding_df %>% dplyr::left_join(disease_transform,by="Tumor_Name")
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<-subset_embedding_df[,c(-1)]
    
    result_run_aff<-run_affinity_clustering(Clusters_Names,subset_embedding_df)
    
    flag_converge <- result_run_aff[[1]]
    subset_affinity_df<-result_run_aff[[2]]
    
    if(flag_converge=="No"){
      for (iter_nested_affinity_cluser in 1: dim(subset_affinity_df)[1]){
        ind_location <- which (affinity_cluster_df$Tumor_Names==subset_affinity_df$Tumor_Names[iter_nested_affinity_cluser])
        affinity_cluster_df$Cluster_ID[ind_location]<-subset_affinity_df$SubCluster_ID[iter_nested_affinity_cluser]
      }
    }else if(flag_converge=="Yes"){
      converge_list<-append(Clusters_Names,converge_list)
    }
    
    
  }
  cluster_frequency_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
  colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
  cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
  large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]
  large_cluster_labels<-setdiff(large_cluster_labels, unlist(converge_list))
}


save(affinity_cluster_df,file = paste(intermediate_dir,"/affinity_cluster_df_ada2.RData",sep=""))

#save.image(file = "script6a.RData")



