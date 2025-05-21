# Computes Kmeans cluster of ADA2 data and also computes silhouette index
start_time<-Sys.time()
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  # library(qdapRegex)
  library(jsonlite)
  library(httr)
  #  library(biomaRt)
  #  library(ghql)
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
plots_dir<-file.path(root_dir,"plots_5th")


# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2_5thed.csv",sep="") )
colnames(disease_transform)[1]<-"Diseases"
rownames(disease_transform)<-disease_transform$Diseases # Needed for AP Clust


# Set Seed
set.seed(13)

ncol= dim(disease_transform)[2]
silhouette_score <- function(k){
  print(k)
  km <- kmeans(disease_transform[,2:ncol], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform[,2:ncol]))
  return(mean(ss[, 3]))
}

k <- c(10,100,500,1000,2000,3000,4000,5000,
       5500,5900, 6000,6050,6100,6200,6500,6800,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)



k1<-c(10,100,500,1000,2000,3000,4000,5000)
k2<-c( 5500,5900, 6000,6050,6100,6200,6500,6800,7000)
k3<-c(8000,9000,10000,11000,12000,13000,14000,15000,16000)


end_time<-Sys.time()

start_time_1<-Sys.time()
avg_sil_1 <- sapply(k1, silhouette_score)#7:37 pm
end_time_1<-Sys.time()
save.image(file = "5a_5th.RData")

start_time_2<-Sys.time()
avg_sil_2 <- sapply(k2, silhouette_score)#7:37 pm
end_time_2<-Sys.time()
save.image(file = "5a_5th.RData")

start_time_3<-Sys.time()
avg_sil_3 <- sapply(k3, silhouette_score)#7:37 pm
end_time_3<-Sys.time()
save.image(file = "5a_5th.RData")

avg_sil<-c(avg_sil_1,avg_sil_2,avg_sil_3)

# avg_sil <- sapply(k, silhouette_score)#7:37 pm
# 
start_time_4<-Sys.time()
Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #6050
# 
Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]
# 
# end_time_5a_5th_1<-Sys.time()
# elapsed_time_5a_5th_1<- as.numeric(difftime(end_time_5a_5th_1, start_time_5a_5th_1, units = "secs"))
# save.image(file = "5a_5th.RData")
# 
# 
# # Kmeans optimal cluster is 6050
# start_time_5a_5th_2<-Sys.time()
# 
km.res <- eclust(disease_transform[,2:ncol], "kmeans", k = Kmeans_silhouette_Max$k,nstart = 25, graph = FALSE)
kmeans_clust_result <- as.data.frame(km.res$cluster)
kmeans_clust_result$Tumors<-rownames(kmeans_clust_result)
colnames(kmeans_clust_result)[1]<-"cluster"
kmeans_clust_result <- kmeans_clust_result %>% dplyr::select(Tumors,cluster)
# 
sil <- silhouette(km.res$cluster, dist(disease_transform[,2:ncol])) # Verify this 
sil<-as.data.frame(sil)
sil$Tumors<-names(km.res$cluster)
# 
kmeans_clust_result <- kmeans_clust_result %>% dplyr::left_join(sil,by=c("cluster", "Tumors"))
# 
kmeans_clust_result<-kmeans_clust_result[order(kmeans_clust_result$cluster),]
rownames(kmeans_clust_result)<-NULL
# 
mean_freq_kmeans <- kmeans_clust_result %>% dplyr::select(cluster, sil_width) %>% dplyr::group_by(cluster) %>% dplyr::summarise(mean_silo_score=mean(sil_width),cluster_member_count =dplyr::n()) 
kmeans_clust_result<- kmeans_clust_result %>% dplyr::left_join(mean_freq_kmeans,by="cluster")
# 
# 
# 
# 
# 
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
ct_tumor_df<-ct_tumor_df[,c(1,2)]
# 
colnames(kmeans_clust_result)[1]<-"Tumor_Names"
colnames(ct_tumor_df)[2]<-"Tumor_Names"
# 
kmeans_clust_result<-kmeans_clust_result%>%left_join(ct_tumor_df,by="Tumor_Names")
# 
kmeans_clust_result<-kmeans_clust_result[,c(7,1:6)]
end_time_4<-Sys.time()

elapsed_time_1<- as.numeric(difftime(end_time_1, start_time_1, units = "secs"))
elapsed_time_2<- as.numeric(difftime(end_time_2, start_time_2, units = "secs"))
elapsed_time_3<- as.numeric(difftime(end_time_3, start_time_3, units = "secs"))
elapsed_time_4<- as.numeric(difftime(end_time_4, start_time_4, units = "secs"))

total_time_elapsed<- elapsed_time_1+elapsed_time_2+elapsed_time_3+elapsed_time_4

objs <- ls(envir = .GlobalEnv)

# Get the size of each object
sizes <- sapply(objs, function(x) object.size(get(x, envir = .GlobalEnv)))

# Convert sizes to readable format and summarize
sizes_df <- data.frame(
  Object = objs,
  Size_MB = round(sizes / (1024^2), 3)
)

# Sort by size descending
sizes_df <- sizes_df[order(-sizes_df$Size_MB), ]

# Print the total memory used
cat("Total memory used in Global Environment:", round(sum(sizes) / (1024^2), 3), "MB\n")

Runtime_5A_5th= round(total_time_elapsed, 2)
Memory_gb_5A_5th=round(sum(sizes) / (1024^2), 3)/1000
# 
# write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding_ada2_5thed.csv",sep=""))
# 
# 
# end_time_5a_5th_2<-Sys.time()#2025-May-19 00:53:59
# 
# elapsed_time_5a_5th_2<- as.numeric(difftime(end_time_5a_5th_2, start_time_5a_5th_2, units = "secs"))
# 
# write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding_ada2.csv",sep=""))
# 
# 
# objs <- ls(envir = .GlobalEnv)
# 
# # Get the size of each object
# sizes <- sapply(objs, function(x) object.size(get(x, envir = .GlobalEnv)))
# 
# # Convert sizes to readable format and summarize
# sizes_df <- data.frame(
#   Object = objs,
#   Size_MB = round(sizes / (1024^2), 3)
# )
# 
# # Sort by size descending
# sizes_df <- sizes_df[order(-sizes_df$Size_MB), ]
# 
# # Print the total memory used
# cat("Total memory used in Global Environment:", round(sum(sizes) / (1024^2), 3), "MB\n")
# 
# time_elapsed_5th<-elapsed_time_5a_5th_1+elapsed_time_5a_5th_2
# Runtime_5A_5th = round(elapsed_time_5a_5th, 2)
# Memory_gb_5a_5th=round(sum(sizes) / (1024^2), 3)/1000
# 
# save.image(file = "5a_5th.RData")
# 
save.image(paste(analysis_dir,"/run-time-analysis/Kmeans_ADA2/5a_5th.RData",sep=""))
