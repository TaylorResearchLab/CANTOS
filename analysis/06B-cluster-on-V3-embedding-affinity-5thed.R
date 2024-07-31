# Computes affinity cluster of V3 data. Nested clustering is performed on large cluster. Cluster size is determined to be large using Z scores on cluster membership.

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
intermediate_dir <- file.path(analysis_dir,"intermediate_5th")
results_dir <- file.path(analysis_dir,"results_5th")

source(paste(util_dir,"/run_affinity_clustering.R",sep=""))


disease_transform_v3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3_5thed.csv",sep="") )
colnames(disease_transform_v3)[1]<-"Tumor_Names"
rownames(disease_transform_v3)<-disease_transform_v3$Tumor_Name # Needed for AP Clust


# Set Seed
set.seed(13)
#affinity cluster


#######V3

dist_euclidean_v3<- dist(disease_transform_v3,method = "euclidean")
dist_euclidean_v3<-as.matrix(dist_euclidean_v3)
simmilarity_euclidean_v3<- 1/(1+dist_euclidean_v3)
af_clust_euclidean_v3 <- apcluster(simmilarity_euclidean_v3)#10:17 pm

