# Computes Kmeans cluster of ADA2 data and also computes silhouette index
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
  library(factoextra)
  library(cluster)
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
result_dir <-file.path(analysis_dir,"results_5th")
plots_dir<-file.path(root_dir,"plots")


# Load PCA Embeddings of CT , WHO, NCIT
disease_transform_v3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3_5thed.csv",sep="") )
colnames(disease_transform_v3)[1]<-"Diseases"
rownames(disease_transform_v3)<-disease_transform_v3$Diseases # Needed for AP Clust

ncol= dim(disease_transform_v3)[2]


# Set Seed
set.seed(13)


# Peform Clustering 

silhouette_score <- function(k){
  km <- kmeans(disease_transform_v3[,2:ncol], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform_v3[,2:ncol]))
  print(k)
  return(mean(ss[, 3]))
}

k <- c(10,500,1000,2000,3000,4000,5000,5500,5800,6000,6100,6200,6300,6400,
       6500,6600,6700,6800,6900,
       7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)

avg_sil <- sapply(k, silhouette_score)
