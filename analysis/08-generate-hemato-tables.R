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
  library(doParallel)
  library(foreach)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))


# 
#load(paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))
affinity_cluster_V3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))



###### Lymphoma Luek analysis

lymphoma_leukemia_strings <- c("leukemia", "lymphoma", "leukemias", "lymphomas", "leukaemia", "leukaemias",
                               "leuk","hematologic tumors","hemato")

affinity_cluster_hema_ADA2_df <- affinity_cluster_ADA2_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))
hema_label <- unique(affinity_cluster_hema_ADA2_df$Cluster_ID)
affinity_cluster_hema_ADA2_df <- affinity_cluster_ADA2_df %>% dplyr::filter(Cluster_ID %in% hema_label )


affinity_cluster_hema_V3_df <- affinity_cluster_V3_df %>% dplyr::filter(str_detect(Tumor_Names,paste(lymphoma_leukemia_strings, collapse = "|")))
hema_label <- unique(affinity_cluster_hema_V3_df$Cluster_ID)
affinity_cluster_hema_V3_df <- affinity_cluster_V3_df %>% dplyr::filter(Cluster_ID %in% hema_label )


#### ISOLATION
affinity_cluster_outlier<-affinity_cluster_hema_ADA2_df%>%dplyr::select(Tumor_Names,Cluster_ID)

disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Tumor_Names"
rownames(disease_transform)<-disease_transform$Tumor_Names 


hemato_embedding<-affinity_cluster_outlier %>% dplyr::left_join(disease_transform,by="Tumor_Names")

affinity_cluster_outlier$isolation_outlier_score<-NA

hemato_cluster_labels <- unique(hemato_embedding$Cluster_ID)
for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
  model <- isolation.forest(hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], ndim=3, ntrees=50, nthreads=1)
  scores <- predict(model, hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], type="score")
  ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
  affinity_cluster_outlier$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
    affinity_cluster_outlier$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_outlier<- affinity_cluster_outlier %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


affinity_cluster_outlier$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  ind_clust <- which(affinity_cluster_outlier$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(hemato_embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(hemato_embedding_subset[,3:137],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_outlier$LOF_Scores[ind_clust]<-0
  }
  

}

affinity_cluster_outlier<- affinity_cluster_outlier %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))

affinity_cluster_hema_df2<- affinity_cluster_hema_df %>% dplyr::select(Tumor_Names,assigned_class,Cluster_ID)
#affinity_cluster_hema_df2 <- affinity_cluster_hema_df2 %>% filter(!Tumor_Names %in% unique(c(NCIT_embedding_df$Disease,WHO_embedding_df$Disease)))
affinity_cluster_hema_df2 <- affinity_cluster_hema_df2 %>% dplyr::left_join(affinity_cluster_outlier[,c(1)],by="Tumor_Names")





##### ISOLATION FOR HEMA DF2
affinity_cluster_outlier2<-affinity_cluster_hema_df2%>%dplyr::select(Tumor_Names,assigned_class)




hemato_embedding<-affinity_cluster_outlier2 %>% dplyr::left_join(disease_transform,by="Tumor_Names")

affinity_cluster_outlier2$isolation_outlier_score<-NA

hemato_cluster_labels <- unique(hemato_embedding$assigned_class)
for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(assigned_class==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], ndim=3, ntrees=50, nthreads=1)
    scores <- predict(model, hemato_embedding_subset[1:nrow(hemato_embedding_subset),3:ncol(hemato_embedding_subset)], type="score")
    ind_clust <- which(affinity_cluster_hema_df2$assigned_class==cluster_label_current)
    affinity_cluster_outlier2$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_hema_df2$assigned_class==cluster_label_current)
    affinity_cluster_outlier2$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_outlier2<- affinity_cluster_outlier2 %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


affinity_cluster_outlier2$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  ind_clust <- which(affinity_cluster_outlier2$assigned_class==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  set.seed(13)
  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(assigned_class==cluster_label_current)
  if(dim(hemato_embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(hemato_embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(hemato_embedding_subset[,3:137],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    affinity_cluster_outlier2$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    affinity_cluster_outlier2$LOF_Scores[ind_clust]<-0
  }
  
  
}


affinity_cluster_outlier2<- affinity_cluster_outlier2 %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))

ind_both <- which(affinity_cluster_outlier2$assigned_class=="Both")
for(iter in ind_both){
  affinity_cluster_outlier2$assigned_class[ind_both]<-affinity_cluster_outlier2$Tumor_Names[ind_both]
}

ind_both <- which(hemato_embedding$assigned_class=="Both")
for(iter in ind_both){
  hemato_embedding$assigned_class[ind_both]<-hemato_embedding$Tumor_Names[ind_both]
}



hemato_cluster_labels<-unique(affinity_cluster_outlier2$assigned_class)
affinity_cluster_outlier2$Z_score<-NA

for(iter in 1:length(hemato_cluster_labels)){
  cluster_label_current <- hemato_cluster_labels[iter]
  ind_clust <- which(affinity_cluster_outlier2$assigned_class==cluster_label_current)

  hemato_embedding_subset <- hemato_embedding %>% dplyr::filter(assigned_class==cluster_label_current)
  
  if(dim(hemato_embedding_subset)[1]>1){
  current_cluster_label_embedding<- current_cluster_label_embedding<- disease_transform[which(rownames(disease_transform)==cluster_label_current),]
  current_cluster_label_embedding<-current_cluster_label_embedding[,c(-1)]
  distance_matrix <- rbind(current_cluster_label_embedding,hemato_embedding_subset[,3:ncol(hemato_embedding_subset)])
  rownames(distance_matrix)<-NULL
  distance_matrix <- as.matrix(dist(distance_matrix,method = "euclidean",diag = TRUE,upper = TRUE))
  distance_matrix<-distance_matrix[c(-1),c(1)]
  distance_matrix<- (distance_matrix-mean(distance_matrix))/sd(distance_matrix)
  affinity_cluster_outlier2$Z_score[ind_clust]<-distance_matrix
  }
  else{
    affinity_cluster_outlier2$Z_score[ind_clust]<-0
    
  }
  }
  
  
affinity_cluster_outlier2<-affinity_cluster_outlier2 %>% dplyr::mutate(Z_score_outlier =case_when(Z_score>2.5 ~ "Yes", TRUE~
                                                                                                    "No"))
affinity_cluster_outlier3<-affinity_cluster_outlier2 %>% dplyr::select(Tumor_Names,assigned_class,Z_score_outlier,Isolation_Outlier,LOF_Outlier)

write.csv(affinity_cluster_hema_df,"hemato_tumor.csv")

stopCluster(cl)
save.image("script08.RData")