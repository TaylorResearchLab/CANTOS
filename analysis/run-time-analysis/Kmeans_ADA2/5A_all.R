# Computes Kmeans cluster of ADA2 data and also computes silhouette index
start_time_5a_all<-Sys.time()
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
#root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
root_dir <- "/home/lahiria/CANTOS-RUN-TIME"
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
result_dir <-file.path(analysis_dir,"results")
plots_dir<-file.path(root_dir,"plots")

print("Loading data")

# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca_ada2.csv",sep="") )
colnames(disease_transform)[1]<-"Diseases"
rownames(disease_transform)<-disease_transform$Diseases # Needed for AP Clust




print("Loading Complete")

# Set Seed
set.seed(13)


# Peform Clustering 
ncol= dim(disease_transform)[2]
silhouette_score <- function(k){
  print(k)
  km <- kmeans(disease_transform[,2:ncol], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform[,2:ncol]))
  return(mean(ss[, 3]))
}

k <- c(10,100,500,1000,2000,3000,4000,5000,
       5500,5900, 6000,6050,6100,6200,6500,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)
print("calculate sil scores ")
avg_sil <- sapply(k, silhouette_score)#2:23pm-5:12

Kmeans_silhouette<-as.data.frame(cbind(k,avg_sil))
colnames(Kmeans_silhouette) <- c("k","mean_silhouette_score") #6200

Kmeans_silhouette_Max <- Kmeans_silhouette[ which(max(Kmeans_silhouette$mean_silhouette_score) == Kmeans_silhouette$mean_silhouette_score), ]



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


ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
ct_tumor_df<-ct_tumor_df[,c(1,2)]

colnames(kmeans_clust_result)[1]<-"Tumor_Names"
colnames(ct_tumor_df)[2]<-"Tumor_Names"



kmeans_clust_result<-kmeans_clust_result%>%left_join(ct_tumor_df,by="Tumor_Names")

kmeans_clust_result<-kmeans_clust_result[,c(7,1:6)]

end_time_5a_all<-Sys.time()

write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding_ada2.csv",sep=""))

elapsed_time_5a_all<- as.numeric(difftime(end_time_5a_all, start_time_5a_all, units = "secs"))

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

Runtime_5A_all = round(elapsed_time_5a_all, 2)
Memory_gb_5a_all=round(sum(sizes) / (1024^2), 3)/1000

save.image(file = "5a_all.RData")
