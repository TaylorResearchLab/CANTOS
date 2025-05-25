#Detect if Affinity cluster members are outliers using LOF and Isolation Forest on ADA2 and V3 data
start_time_1<-Sys.time()

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
print("Running 8_all")

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")
intermediate_dir <- file.path(run_time_analysis,"intermediate")


# Load affinity clusters for ADA2 
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep="")) # 10 cols
affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df[,c(-1)]

# Load embeddings 
disease_transform_ADA2<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2.csv",sep="") )
colnames(disease_transform_ADA2)[1]<-"Tumor_Names"
rownames(disease_transform_ADA2)<-disease_transform_ADA2$Tumor_Names 

print("CPT1")


# Compute isolation forest
set.seed(13)

embedding_ADA2<-affinity_cluster_ADA2_df %>% dplyr::left_join(disease_transform_ADA2,by="Tumor_Names")
affinity_cluster_ADA2_df$isolation_outlier_score<-NA
idx_ada2<-which(colnames(embedding_ADA2)=="PC1")

cluster_labels_ADA2 <- unique(embedding_ADA2$Cluster_ID)
for(iter in 1:length(cluster_labels_ADA2)){
  cluster_label_current <- cluster_labels_ADA2[iter]
  embedding_subset <- embedding_ADA2 %>% dplyr::filter(Cluster_ID==cluster_label_current)
  if(dim(embedding_subset)[1]>2){ # Need at least 2 data points to run isolation forest
    model <- isolation.forest(embedding_subset[1:nrow(embedding_subset),idx_ada2:ncol(embedding_subset)], ndim=3, ntrees=100, nthreads=1) # ntrees 50 initially
    scores <- predict(model, embedding_subset[1:nrow(embedding_subset),idx_ada2:ncol(embedding_subset)], type="score")
    ind_clust <- which(affinity_cluster_ADA2_df$Cluster_ID==cluster_label_current)
    affinity_cluster_ADA2_df$isolation_outlier_score[ind_clust]<-scores
  }else{
    ind_clust <- which(affinity_cluster_ADA2_df$Cluster_ID==cluster_label_current)
    affinity_cluster_ADA2_df$isolation_outlier_score[ind_clust]<-0
  }
}
affinity_cluster_ADA2_df<- affinity_cluster_ADA2_df %>% dplyr::mutate(Isolation_Outlier = case_when(isolation_outlier_score>0.5 ~ "Yes", TRUE ~ "No"))
print("CPT2")

stop_time_1<-Sys.time()
save.image(paste(run_time_analysis,"/ADA002_AP/8-ADA002_all.RData",sep=""))

start_time_2<-Sys.time()

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
      lof_scores_minpts <- lof(embedding_subset[,idx_ada2:ncol(embedding_subset)],iter_pts)
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

print("CPT3")


stop_time_2<-Sys.time()

elapsed_time_1<- as.numeric(difftime(stop_time_1, start_time_1, units = "secs"))
elapsed_time_2<- as.numeric(difftime(stop_time_2, start_time_2, units = "secs"))

elapsed_time<-sum(elapsed_time_1,elapsed_time_2)




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

Runtime_8_all = round(elapsed_time, 2)
Memory_gb_8_all=round(sum(sizes) / (1024^2), 3)/1000



####
write.csv(affinity_cluster_ADA2_df,paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))


save.image(paste(run_time_analysis,"/ADA002_AP/8-ADA002_all.RData",sep=""))
