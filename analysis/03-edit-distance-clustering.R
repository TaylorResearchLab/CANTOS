# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(apcluster)
  library(stringdist)
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

# Load data

load(paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


cluster_results_lv<- read.csv(paste(results_dir,"/cluster_lv.csv",sep=""))
cluster_results_jw<- read.csv(paste(results_dir,"/cluster_jw.csv",sep=""))
cluster_results_cosine<- read.csv(paste(results_dir,"/cluster_cosine.csv",sep=""))


# Kmeans 

simmilarity_matrix_cosine = 1 - dissimilarity_matrix_cosine
simmilarity_matrix_jw = 1-dissimilarity_matrix_jw


apres1b <- apcluster(1-dissimilarity_matrix_jw[1:100,1:100])
affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apres1b@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(apres1b@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
