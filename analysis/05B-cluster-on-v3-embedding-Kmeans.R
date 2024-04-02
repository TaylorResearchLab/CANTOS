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
intermediate_dir <- file.path(analysis_dir,"intermediate")
result_dir <-file.path(analysis_dir,"results")
plots_dir<-file.path(root_dir,"plots")


# Load PCA Embeddings of CT , WHO, NCIT
disease_transform_v3<- read.csv(paste(intermediate_dir,"/disease_transform_pca_v3.csv",sep="") )
colnames(disease_transform_v3)[1]<-"Diseases"
rownames(disease_transform_v3)<-disease_transform_v3$Diseases # Needed for AP Clust



# Set Seed
set.seed(13)


# Peform Clustering 

silhouette_score <- function(k){
  km <- kmeans(disease_transform_v3[,2:178], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform_v3[,2:136]))
  return(mean(ss[, 3]))
}

k <- c(10,500,1000,3000,5000,6000,8000,10000,13000,15000,16000)
avg_sil <- sapply(k, silhouette_score)#11:04-12:18

Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #6000

Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]

