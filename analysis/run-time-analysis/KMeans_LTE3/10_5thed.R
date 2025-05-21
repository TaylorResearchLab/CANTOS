start_time<-Sys.time()
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(readxl)
  library(dbscan)
  library(isotree)
  library(qdapRegex)
  library(ghql)
})

setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")
intermediate_dir <- file.path(run_time_analysis,"intermediate_5th")
result_dir <-file.path(run_time_analysis,"results_5th")

source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))



#Read Kmeans 
kmeans_clust_result_embedding_V3 <- read_csv(paste(result_dir,"/kmeans_clust_result_embedding_v3_5thed.csv",sep=""))
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3[,c(-1)]
who_ncit_match_v3<-read.csv(paste(analysis_dir,"/results_5th/who_ncit_match_v3_5thed.csv",sep=""))
who_ncit_match_v3<-who_ncit_match_v3[,c(-1)]


# Join the matches to Kmeans 
kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::select(nct_id,Tumor_Names,cluster)

colnames(kmeans_clust_result_embedding_V3)[3]<-c("Cluster_ID")

kmeans_clust_result_embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::left_join(who_ncit_match_v3,by="Tumor_Names")

kmeans_clust_result_embedding_V3<- cluster_label_assignment_refined(kmeans_clust_result_embedding_V3)


disease_transform_V3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3_5thed.csv",sep="") )
colnames(disease_transform_V3)[1]<-"Tumor_Names"
rownames(disease_transform_V3)<-disease_transform_V3$Tumor_Names 


set.seed(13)

### Compute isolation forest for V3 embedding Kmeans

# Compute isolation forest
embedding_V3<-kmeans_clust_result_embedding_V3 %>% dplyr::left_join(disease_transform_V3,by="Tumor_Names")
kmeans_clust_result_embedding_V3$isolation_outlier_score<-NA
idx_v3<-which(colnames(embedding_V3)=="PC1")

cluster_labels_V3 <- unique(embedding_V3$Cluster_ID)
for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),idx_v3:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),idx_v3:ncol(embedding_subset)], type="score")
    ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_V3$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
    kmeans_clust_result_embedding_V3$isolation_outlier_score[ind_clust]<-0
  }
}
kmeans_clust_result_embedding_V3<- kmeans_clust_result_embedding_V3 %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))


# Compute LOF 
kmeans_clust_result_embedding_V3$LOF_Scores<-NA
lof_scores_minpts_list<-list()

for(iter in 1:length(cluster_labels_V3)){
  cluster_label_current <- cluster_labels_V3[iter]
  ind_clust <- which(kmeans_clust_result_embedding_V3$Cluster_ID==cluster_label_current)
  lof_scores_minpts_list<-list()
  
  embedding_subset <- embedding_V3 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    min_pts<- 2:(dim(embedding_subset)[1]-1)
    for(iter_pts in min_pts){
      lof_scores_minpts <- lof(embedding_subset[,idx_v3:ncol(embedding_subset)],iter_pts)
      lof_scores_minpts_list[[as.character(iter_pts)]]<-lof_scores_minpts
    }
    lof_scores_minpts_list<- t(as.data.frame(lof_scores_minpts_list))
    lof_scores_minpts_list_median<-apply(lof_scores_minpts_list,2,median)
    kmeans_clust_result_embedding_V3$LOF_Scores[ind_clust]<-lof_scores_minpts_list_median
    
  }else{
    kmeans_clust_result_embedding_V3$LOF_Scores[ind_clust]<-0
  }
  
  
}

kmeans_clust_result_embedding_V3<- kmeans_clust_result_embedding_V3 %>% dplyr::mutate(LOF_Outlier = case_when(LOF_Scores>1 ~ "Yes", TRUE ~ "No"))
kmeans_clust_result_embedding_V3_short <- kmeans_clust_result_embedding_V3 %>% dplyr::select(nct_id,Tumor_Names,who_cluster_label)
colnames(kmeans_clust_result_embedding_V3_short)[3]<-"kmeans_v3"
kmeans_clust_result_embedding_V3_short_NCIT <- kmeans_clust_result_embedding_V3 %>% dplyr::select(nct_id,Tumor_Names,ncit_cluster_label)
colnames(kmeans_clust_result_embedding_V3_short_NCIT)[3]<-"kmeans_v3"

end_time<-Sys.time()

elapsed_time<- as.numeric(difftime(end_time, start_time, units = "secs"))
objs <- ls(envir = .GlobalEnv)

# Get the size of each object
sizes <- sapply(objs, function(x) object.size(get(x, envir = .GlobalEnv)))

# Convert sizes to readable format and summarize
sizes_df <- data.frame(
  Object = objs,
  Size_MB = round(sizes / (1024^2), 3)
)

# Sort by size descending
sizes_df <- sizes_df[order(-sizes_df$Size_MB), ]

# Print the total memory used
cat("Total memory used in Global Environment:", round(sum(sizes) / (1024^2), 3), "MB\n")

save.image(paste(run_time_analysis,"/KMeans_LTE3/10_5th.RData",sep=""))


