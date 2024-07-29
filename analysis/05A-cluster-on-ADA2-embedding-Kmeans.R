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
intermediate_dir <- file.path(analysis_dir,"intermediate")
result_dir <-file.path(analysis_dir,"results")
plots_dir<-file.path(root_dir,"plots")


# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2.csv",sep="") )
colnames(disease_transform)[1]<-"Diseases"
rownames(disease_transform)<-disease_transform$Diseases # Needed for AP Clust


# Set Seed
set.seed(13)


# Peform Clustering 
ncol= dim(disease_transform)[2]
silhouette_score <- function(k){
  km <- kmeans(disease_transform[,2:ncol], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform[,2:ncol]))
  return(mean(ss[, 3]))
}

k <- c(10,100,500,1000,2000,3000,4000,5000,
       5500,5900, 6000,6050,6100,6200,6500,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)
avg_sil <- sapply(k, silhouette_score)#2:23pm-5:12

Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #6200

Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]


p1<-ggplot(Kmeans_silhouette, aes(x =k, y = mean_silhouette_score)) + geom_point() +
  geom_point(data = Kmeans_silhouette[which.max(Kmeans_silhouette$mean_silhouette_score), ], color="red")+
  scale_x_continuous("k", labels = as.character(k), breaks = k) + ggtitle("Kmean Silhouette Score vs Clusters")


# Kmeans optimal cluster is 6000
km.res <- eclust(disease_transform[,2:ncol], "kmeans", k = Kmeans_silhouette_Max$k,nstart = 25, graph = FALSE)
kmeans_clust_result <- as.data.frame(km.res$cluster)
kmeans_clust_result$Tumors<-rownames(kmeans_clust_result)
colnames(kmeans_clust_result)[1]<-"cluster"
kmeans_clust_result <- kmeans_clust_result %>% dplyr::select(Tumors,cluster)

sil <- silhouette(km.res$cluster, dist(disease_transform[,2:ncol])) # Verify this 
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

ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
ct_tumor_df<-ct_tumor_df[,c(1,2)]

colnames(kmeans_clust_result)[1]<-"Tumor_Names"
colnames(display_table_benchmark_kmeans)[1]<-"Tumor_Names"
colnames(ct_tumor_df)[2]<-"Tumor_Names"

kmeans_clust_result<-kmeans_clust_result%>%left_join(ct_tumor_df,by="Tumor_Names")
display_table_benchmark_kmeans<-display_table_benchmark_kmeans%>%left_join(ct_tumor_df,by="Tumor_Names")

kmeans_clust_result<-kmeans_clust_result[,c(7,1:6)]
display_table_benchmark_kmeans<-display_table_benchmark_kmeans[,c(7,1:6)]

write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding.csv",sep=""))
write.csv(display_table_benchmark_kmeans,paste(result_dir,"/display_table_benchmark_kmeans.csv",sep=""))


# Plot for Kmeans
clust_plot_kmeans <- display_table_benchmark_kmeans %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) 
clust_plot_kmeans<-unique(clust_plot_kmeans)
p_kmeans_benchmark <- ggplot(clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.005) + labs(title =  paste("Kmeans clusters, K=",Kmeans_silhouette_Max$k,sep=""))
ggsave(p_kmeans_benchmark, filename = paste(plots_dir,"/kmeans_Embedding_Benchmark_ada2.png",sep=""), height = 30, width = 21, units = "cm")

# Global Plot kmeans
global_clust_plot_kmeans <- kmeans_clust_result %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
plt_global_kmeans <- ggplot(global_clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title =  paste("Kmeans clusters, K=",Kmeans_silhouette_Max$k,sep=""))
ggsave(plt_global_kmeans, filename = paste(plots_dir,"/plt_global_kmeans_embedding_ada2.pdf",sep=""), height = 30, width = 21, units = "cm")



save.image(file = "script5a.RData")
