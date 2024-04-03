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

k <- c(10,100,500,1000,2000,3000,4000,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,
       5500,5900, 6000,6050,6100,6200,6500,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)

avg_sil <- sapply(k, silhouette_score)

Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #5400

Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]

p1<-ggplot(Kmeans_silhouette, aes(x =k, y = mean_silhouette_score)) + geom_point() +
  geom_point(data = Kmeans_silhouette[which.max(Kmeans_silhouette$mean_silhouette_score), ], color="red")+
  scale_x_continuous("k", labels = as.character(k), breaks = k) + ggtitle("Kmean Silhouette Score vs Clusters with V3 embeddings")


# Kmeans optimal cluster is 5000
km.res <- eclust(disease_transform_v3[,2:ncol(disease_transform_v3)], "kmeans", k = Kmeans_silhouette_Max$k,nstart = 25, graph = FALSE)
kmeans_clust_result <- as.data.frame(km.res$cluster)
kmeans_clust_result$Tumors<-rownames(kmeans_clust_result)
colnames(kmeans_clust_result)[1]<-"cluster"
kmeans_clust_result <- kmeans_clust_result %>% dplyr::select(Tumors,cluster)

sil <- silhouette(km.res$cluster, dist(disease_transform_v3[,2:ncol(disease_transform_v3)])) # Verify this 
sil<-as.data.frame(sil)
sil$Tumors<-names(km.res$cluster)

kmeans_clust_result <- kmeans_clust_result %>% dplyr::left_join(sil,by=c("cluster", "Tumors"))

kmeans_clust_result<-kmeans_clust_result[order(kmeans_clust_result$cluster),]
rownames(kmeans_clust_result)<-NULL

mean_freq_kmeans <- kmeans_clust_result %>% dplyr::select(cluster, sil_width) %>% dplyr::group_by(cluster) %>% dplyr::summarise(mean_silo_score=mean(sil_width),cluster_member_count =dplyr::n()) 
kmeans_clust_result<- kmeans_clust_result %>% dplyr::left_join(mean_freq_kmeans,by="cluster")


# Show results with following tumors 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

cluster_ind_benchmark_tumor <- kmeans_clust_result$cluster[kmeans_clust_result$Tumors %in% benchmark_tumors]

display_table_benchmark_kmeans <- kmeans_clust_result %>% filter(cluster %in% cluster_ind_benchmark_tumor)
display_table_benchmark_kmeans<- display_table_benchmark_kmeans[order(display_table_benchmark_kmeans$cluster),]
rownames(display_table_benchmark_kmeans)<-NULL



write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding_v3.csv",sep=""))
write.csv(display_table_benchmark_kmeans,paste(result_dir,"/display_table_benchmark_kmeans_v3.csv",sep=""))


# Plot for Kmeans
clust_plot_kmeans <- display_table_benchmark_kmeans %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) 
clust_plot_kmeans<-unique(clust_plot_kmeans)
p_kmeans_benchmark <- ggplot(clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.005) + labs(title = "Kmeans clusters, K=5400")
ggsave(p_kmeans_benchmark, filename = paste(plots_dir,"/kmeans_Embedding_Benchmark_v3.png",sep=""), height = 30, width = 21, units = "cm")

# Global Plot kmeans
global_clust_plot_kmeans <- kmeans_clust_result %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
plt_global_kmeans <- ggplot(global_clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Kmeans clusters, K=5400")
ggsave(plt_global_kmeans, filename = paste(plots_dir,"/plt_global_kmeans_embedding_v3.pdf",sep=""), height = 30, width = 21, units = "cm")

