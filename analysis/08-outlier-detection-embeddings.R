#Detect if Affinity cluster members are outliers using LOF and Isolation Forest on ADA2 and V3 data
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(jsonlite)
  library(ghql)
  library(readxl)
  library(doParallel)
  library(foreach)
  library(isotree)
  library(dbscan)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")


# Load affinity clusters for ADA2 and V3
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep="")) # 8 cols
affinity_cluster_v3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep="")) # 8 cols

affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df[,c(-1)]
affinity_cluster_v3_df<- affinity_cluster_v3_df[,c(-1)]

# Load embeddings 
disease_transform_ADA2<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform_ADA2)[1]<-"Tumor_Names"
rownames(disease_transform_ADA2)<-disease_transform_ADA2$Tumor_Names 


disease_transform_V3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
colnames(disease_transform_V3)[1]<-"Tumor_Names"
rownames(disease_transform_V3)<-disease_transform_V3$Tumor_Names 


# Compute isolation forest
set.seed(13)

embedding_ADA2<-affinity_cluster_ADA2_df %>% dplyr::left_join(disease_transform_ADA2,by="Tumor_Names")
affinity_cluster_ADA2_df$isolation_outlier_score<-NA

cluster_labels_ADA2 <- unique(embedding_ADA2$Cluster_ID)
for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], type="score")
    ind_clust <- which(affinity_cluster_ADA2_df$Cluster_ID==cluster_label_current)
    affinity_cluster_ADA2_df$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_ADA2_df$Cluster_ID==cluster_label_current)
    affinity_cluster_ADA2_df$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_ADA2_df<- affinity_cluster_ADA2_df %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


# Compute LOF 
affinity_cluster_ADA2_df$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  ind_clust <- which(affinity_cluster_ADA2_df$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(embedding_subset[,9:ncol(embedding_subset)],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_ADA2_df$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_ADA2_df$LOF_Scores[ind_clust]<-0
  }
  
  
}

affinity_cluster_ADA2_df<- affinity_cluster_ADA2_df %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))



# Compute isolation forest

embedding_V3<-affinity_cluster_v3_df %>% dplyr::left_join(disease_transform_V3,by="Tumor_Names")
affinity_cluster_v3_df$isolation_outlier_score<-NA

cluster_labels_V3 <- unique(embedding_V3$Cluster_ID)
for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),9:ncol(embedding_subset)], type="score")
    ind_clust <- which(affinity_cluster_v3_df$Cluster_ID==cluster_label_current)
    affinity_cluster_v3_df$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_v3_df$Cluster_ID==cluster_label_current)
    affinity_cluster_v3_df$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_v3_df<- affinity_cluster_v3_df %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


# Compute LOF 
affinity_cluster_v3_df$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  ind_clust <- which(affinity_cluster_v3_df$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(embedding_subset[,9:ncol(embedding_subset)],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_v3_df$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_v3_df$LOF_Scores[ind_clust]<-0
  }
  
  
}

affinity_cluster_v3_df<- affinity_cluster_v3_df %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))

####
write.csv(affinity_cluster_ADA2_df,paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))
write.csv(affinity_cluster_v3_df,paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))

