suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
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

# Read file v3
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))
# Read ADA2 
affinity_cluster_ADA2_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_reassigned_df.csv",sep=""))

#Read Kmeans 
kmeans_clust_result_embedding_ADA2 <- read_csv("analysis/results/kmeans_clust_result_embedding.csv")
kmeans_clust_result_embedding_V3 <- read_csv("analysis/results/kmeans_clust_result_embedding_v3.csv")


# Read edit distance cluster
nested_affinity_cluster_cosine <- read_csv("analysis/results/nested_affinity_cluster_cosine.csv")
nested_affinity_cluster_jw <- read_csv("analysis/results/nested_affinity_cluster_jw.csv")
nested_affinity_cluster_lv <- read_csv("analysis/results/nested_affinity_cluster_lv.csv")

