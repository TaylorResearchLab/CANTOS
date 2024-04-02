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

k <- c(10,100,500,1000,2000,3000,4000,5000,
       5500,5900, 6000,6050,6100,6200,6500,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)

avg_sil <- sapply(k, silhouette_score)

Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #6000

Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]

p1<-ggplot(Kmeans_silhouette, aes(x =k, y = mean_silhouette_score)) + geom_point() +
  geom_point(data = Kmeans_silhouette[which.max(Kmeans_silhouette$mean_silhouette_score), ], color="red")+
  scale_x_continuous("k", labels = as.character(k), breaks = k) + ggtitle("Kmean Silhouette Score vs Clusters")


# Kmeans optimal cluster is 6000
km.res <- eclust(disease_transform[,2:136], "kmeans", k = Kmeans_silhouette_Max$k,nstart = 25, graph = FALSE)
kmeans_clust_result <- as.data.frame(km.res$cluster)
kmeans_clust_result$Tumors<-rownames(kmeans_clust_result)
colnames(kmeans_clust_result)[1]<-"cluster"
kmeans_clust_result <- kmeans_clust_result %>% dplyr::select(Tumors,cluster)

sil <- silhouette(km.res$cluster, dist(disease_transform[,2:136])) # Verify this 
sil<-as.data.frame(sil)
sil$Tumors<-names(km.res$cluster)

kmeans_clust_result <- kmeans_clust_result %>% dplyr::left_join(sil,by=c("cluster", "Tumors"))

kmeans_clust_result<-kmeans_clust_result[order(kmeans_clust_result$cluster),]
rownames(kmeans_clust_result)<-NULL

mean_freq_kmeans <- kmeans_clust_result %>% dplyr::select(cluster, sil_width) %>% dplyr::group_by(cluster) %>% dplyr::summarise(mean_silo_score=mean(sil_width),cluster_member_count =dplyr::n()) 
kmeans_clust_result<- kmeans_clust_result %>% dplyr::left_join(mean_freq_kmeans,by="cluster")
