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
# Load affinity Cluster
#load(paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))
source(paste(util_dir,"/nested_affinity_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))



########################################*************************************###########
# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Tumor_Name"
rownames(disease_transform)<-disease_transform$Tumor_Name # Needed for AP Clust


embedding_v3_large <- read.csv(paste(data_dir,"/embedding_tumor_names_text-embedding-3-large_embeddings.csv",sep=""))
colnames(embedding_v3_large)[1]<-"Tumor_Names"
disease_transform_v3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
disease_transform_v3<-disease_transform_v3[,c(-1)]
rownames(disease_transform_v3)<-embedding_v3_large$Tumor_Name # Needed for AP Clust
rm(embedding_v3_large)

# Set Seed
set.seed(13)
#affinity cluster

dist_euclidean<- dist(disease_transform,method = "euclidean")
dist_euclidean<-as.matrix(dist_euclidean)
simmilarity_euclidean<- 1/(1+dist_euclidean)
af_clust_euclidean <- apcluster(simmilarity_euclidean) # 1:24 am - 2:55 am still going....
cat("affinity propogation optimal number of clusters:", length(af_clust_euclidean@clusters), "\n")

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




#################









# Cluster voting
affinity_cluster_nested<- cluster_label_assignment(affinity_cluster_nested)

disease_affinity_cluster_table<- affinity_cluster_nested %>% dplyr::select(Tumor_Names,cluster_label)

# Write
write.csv(affinity_cluster_nested,paste(intermediate_dir,"/affinity_cluster_nested.csv",sep=""))

save.image(file = "script6_affinitycluster.RData")
# write files 
#save(d.apclus2,file = paste(intermediate_dir,"/d.apclus2.RData",sep=""))
save(affinity_cluster_df,file = paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
save(affinity_cluster_annotation,file = paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))

# Silos not computed
source("~/Desktop/MTP_Paper/CT-Embedding-Paper/util/compute_silhouette.R")
affinity_cluster_df2<-affinity_cluster_df
colnames(affinity_cluster_df2)[2]<-"SubsetCluster_IDs"
affinity_cluster_df2<-compute_silhouette(affinity_cluster_df2,dist_euclidean) # Change colname to sublu
save.image(file = "script6_affinitycluster.RData")
save.image(file = "script6_affinitycluster_v3.RData")



mean_freq_af <- affinity_cluster_df2 %>%dplyr::group_by(SubsetCluster_IDs) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
affinity_cluster_df2<- affinity_cluster_df2 %>% dplyr::left_join(mean_freq_af,by="SubsetCluster_IDs")

benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

cluster_ind_benchmark_tumor <- affinity_cluster_df2$SubsetCluster_IDs[affinity_cluster_df2$Tumor_Names %in% benchmark_tumors]

display_table_benchmark_af <- affinity_cluster_df2 %>% filter(SubsetCluster_IDs %in% cluster_ind_benchmark_tumor)
display_table_benchmark_af<- display_table_benchmark_af[order(display_table_benchmark_af$SubsetCluster_IDs),]
rownames(display_table_benchmark_af)<-NULL

write.csv(display_table_benchmark_af,paste(results_dir,"/display_table_benchmark_embedding_af.csv",sep=""))