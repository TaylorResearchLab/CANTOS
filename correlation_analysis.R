#This script is used evaluate the accuracy of each clustering+standardization methods. 

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(ghql)
  library(readxl)
  library(dbscan)
  library(isotree)
  library(infotheo)
  library(combinat)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results")
result_dir_5th <- file.path(analysis_dir,"results_5th")


#tumor_all_gt<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))
#tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))

tumor_all_gt<-read.csv(paste("/Users/lahiria/Desktop/MTP_Paper/CT-Embedding-Paper/analysis/results","/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste("/Users/lahiria/Desktop/MTP_Paper/CT-Embedding-Paper/analysis/results_5th","/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

data_5th<- tumor_5thed_gt[,seq(6,76,2)]
data_all<- tumor_all_gt[,seq(6,76,2)]

col_select <- c("valid_euclidean_dist_v3","valid_euclidean_dist_ada2","valid_af_ad2","valid_af_v3","valid_kmeans_v3",
                "valid_kmeans_ada2","valid_euclidean_dist_MiniLM_L6_v2","valid_euclidean_mpnet_base",
                "valid_euclidean_dist_e5_large","valid_euclidean_dist_MiniLM_L12_v2","valid_euclidean_dist_nomic")


data_5th<-data_5th%>%dplyr::select(col_select)
data_all<-data_all%>%dplyr::select(col_select)

num_methods<-length(col_select)

# Compute pairwise Mutual Information
MI_mat_5th <- matrix(0, nrow = num_methods, ncol = num_methods)
for (i in 1:num_methods) {
  for (j in i:num_methods) {
    mi_score <- mutinformation(data_5th[[i]], data_5th[[j]])
    MI_mat_5th[i, j] <- mi_score
    MI_mat_5th[j, i] <- mi_score  # Fill both halves since it's symmetric
  }
}
colnames(MI_mat_5th) <- colnames(data_5th)
rownames(MI_mat_5th) <- colnames(data_5th)


MI_mat_all <- matrix(0, nrow = num_methods, ncol = num_methods)

# Compute pairwise Mutual Information
for (i in 1:num_methods) {
  for (j in i:num_methods) {
    mi_score <- mutinformation(data_all[[i]], data_all[[j]])
    MI_mat_all[i, j] <- mi_score
    MI_mat_all[j, i] <- mi_score  # Fill both halves since it's symmetric
  }
}
colnames(MI_mat_all) <- colnames(data_all)
rownames(MI_mat_all) <- colnames(data_all)


which(MI_mat_5th==min(MI_mat_5th),arr.ind = TRUE)
which(MI_mat_all==min(MI_mat_all),arr.ind = TRUE)  
method_combinations <- combn(col_select, 3, simplify = FALSE)

compute_avg_mi <- function(methods, mi_matrix) {
  sub_matrix <- mi_matrix[methods, methods]  # Extract sub-matrix for selected methods
  avg_mi <- mean(sub_matrix[lower.tri(sub_matrix)])  # Average off-diagonal mutual information
  return(avg_mi)
}

mi_scores_5th <- sapply(method_combinations, function(combo) compute_avg_mi(combo, MI_mat_5th))
mi_scores_all <- sapply(method_combinations, function(combo) compute_avg_mi(combo, MI_mat_all))

collapsed_string <- sapply(method_combinations, function(x) paste(x, collapse = "+"))


mi_scores_5th<-as.data.frame(cbind(collapsed_string,mi_scores_5th))
mi_scores_all<-as.data.frame(cbind(collapsed_string,mi_scores_all))

colnames(mi_scores_5th)<- c("methods","mi_score")
colnames(mi_scores_all)<- c("methods","mi_score")

mi_scores_5th<-mi_scores_5th[order(mi_scores_5th$mi_score),]
mi_scores_all<-mi_scores_all[order(mi_scores_all$mi_score),]
rownames(mi_scores_5th)<-NULL
rownames(mi_scores_all)<-NULL

print(paste("Best diverse combination (lowest mutual information):", mi_scores_5th$methods[1]))
print(paste("Best diverse combination (lowest mutual information):", mi_scores_all$methods[1]))


