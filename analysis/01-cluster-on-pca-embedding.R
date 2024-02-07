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
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Diseases"


# Peform Clustering 

# Find optimal number of Clusters using KMeans Silhouette 
cluster_results<-fviz_nbclust(disease_transform, kmeans, method = 'silhouette',  k.max = 5000,iter.max=50)
cluster_results_verbose<-fviz_nbclust_verbose(disease_transform, kmeans, method = 'silhouette',  k.max = 13433,iter.max=50)
index_opt_clust<- which(cluster_results$data$y==max(cluster_results$data$y))
opt_clust_size<- as.integer(cluster_results$data$clusters[index_opt_clust]) # 4800
kmeans_disease = kmeans(disease_transform, centers = opt_clust_size, nstart = 100)
diseases_cluster_kmeans <- as.data.frame(kmeans_disease$cluster)


## CHI Index
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")
CH_scroes <- as.data.frame(CH_Results$data$CHIndex) # ratio of the between-cluster variance and the within-cluster variance
WSS_scores<- as.data.frame(CH_Results$data$wss)
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")


#affinity cluster
set.seed(13)
d.apclus2 <- apcluster(negDistMat(r=2), disease_transform)
cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n") #1255

affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(d.apclus2@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(d.apclus2@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
