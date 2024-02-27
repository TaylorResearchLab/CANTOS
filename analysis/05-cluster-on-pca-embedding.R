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
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Diseases"
rownames(disease_transform)<-disease_transform$Diseases # Needed for AP Clust


# Set Seed
set.seed(13)


# Peform Clustering 

silhouette_score <- function(k){
  km <- kmeans(disease_transform[,2:136], centers = k, nstart=25)
  ss <- silhouette(km$cluster, dist(disease_transform[,2:136]))
  return(mean(ss[, 3]))
}

k <- c(10,100,500,1000,2000,3000,4000,5000,
       5500,5900, 6000,6050,6100,6200,6500,
       7000,8000,9000,10000,11000,12000,13000,
       14000,15000,16000)
avg_sil <- sapply(k, silhouette_score)#11:04-12:18

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



# Show results with following tumors 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

cluster_ind_benchmark_tumor <- kmeans_clust_result$cluster[kmeans_clust_result$Tumors %in% benchmark_tumors]

display_table_benchmark_kmeans <- kmeans_clust_result %>% filter(cluster %in% cluster_ind_benchmark_tumor)
display_table_benchmark_kmeans<- display_table_benchmark_kmeans[order(display_table_benchmark_kmeans$cluster),]
rownames(display_table_benchmark_kmeans)<-NULL



write.csv(kmeans_clust_result,paste(result_dir,"/kmeans_clust_result_embedding.csv",sep=""))
write.csv(display_table_benchmark_kmeans,paste(result_dir,"/display_table_benchmark_kmeans.csv",sep=""))


# Plot for Kmeans
clust_plot_kmeans <- display_table_benchmark_kmeans %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) 
clust_plot_kmeans<-unique(clust_plot_kmeans)
p_kmeans_benchmark <- ggplot(clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.005) + labs(title = "Kmeans clusters, K=6000")

# Global Plot kmeans
global_clust_plot_kmeans <- kmeans_clust_result %>% dplyr::select(cluster,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
plt_global_kmeans <- ggplot(global_clust_plot_kmeans, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_kmeans$cluster,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Kmeans clusters, K=6000")
ggsave(plt_global_kmeans, filename = paste(plots_dir,"/plt_global_kmeans_embedding.pdf",sep=""), height = 30, width = 21, units = "cm")





# Find optimal number of Clusters using KMeans Silhouette 
cluster_results<-fviz_nbclust(disease_transform[,2:137], kmeans, method = 'silhouette',  k.max = 5000,iter.max=50)
index_opt_clust<- which(cluster_results$data$y==max(cluster_results$data$y))
opt_clust_size<- as.integer(cluster_results$data$clusters[index_opt_clust]) # 4800
kmeans_disease = kmeans(disease_transform, centers = opt_clust_size, nstart = 100)
diseases_cluster_kmeans <- as.data.frame(kmeans_disease$cluster)
diseases_cluster_kmeans<-cbind(disease_transform$Diseases,diseases_cluster_kmeans)
rownames(diseases_cluster_kmeans)<-NULL






## CHI Index
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")
CH_scroes <- as.data.frame(CH_Results$data$CHIndex) # ratio of the between-cluster variance and the within-cluster variance
WSS_scores<- as.data.frame(CH_Results$data$wss)
CH_Results<-CHCriterion(disease_transform_pca_scaled, kmax=13434,clustermethod="hclust", method = "average")


#affinity cluster
d.apclus2 <- apcluster(negDistMat(r=2), disease_transform) # 1 hr 28 mins
cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n") #1113 clusters 

affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(d.apclus2@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(d.apclus2@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')


# Find cluster membership frequencies 
#affinity_cluster_df$Cluster_Total_Members <- NA
count_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
colnames(count_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
count_table$Cluster_ID<- as.numeric(count_table$Cluster_ID)
affinity_cluster_df <- affinity_cluster_df %>% dplyr::left_join(count_table,by="Cluster_ID")


# Find if tumor has children and pediatric terms in their tumor names and make a seperate column for pediatric cluster
ind_ped <- c(which(str_detect(affinity_cluster_df$Tumor_Names, "childhood")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "children")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "child")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "pediatric")),
             which(str_detect(affinity_cluster_df$Tumor_Names, "paediatric")))

affinity_cluster_df$contains_pediatric_string <- "No"

affinity_cluster_df$contains_pediatric_string[ind_ped]<-"Yes"

pediatric_cluster_ID<-0

affinity_cluster_df$Pediatric_SubsetCluster_ID <- NA  
affinity_cluster_df$Pediatric_SubsetCluster_ID[ind_ped]<-pediatric_cluster_ID


#Nested affinity clusters
barplot(height=affinity_cluster_df$Primary_Cluster_Frequency, names=affinity_cluster_df$Cluster_ID, 
        xlab="cluster_id", 
        ylab="frequency", 
        main="Cluster_Frequency", 
        ylim=c(0,50)
)

large_cluster_index <- which(affinity_cluster_df$Primary_Cluster_Frequency > median(affinity_cluster_df$Primary_Cluster_Frequency))

affinity_cluster_df$High_Primary_Cluster_Membership<-"No"
affinity_cluster_df$High_Primary_Cluster_Membership[large_cluster_index]<-"Yes"

affinity_cluster_df <- affinity_cluster_df %>% dplyr::select(Tumor_Names,Cluster_ID,Primary_Cluster_Frequency,
                                                             High_Primary_Cluster_Membership,contains_pediatric_string,
                                                             Pediatric_SubsetCluster_ID)
Clusters_Names <- unique(affinity_cluster_df$Cluster_ID)

#disease_transform$Tumor_Name<-rownames(disease_transform)
colnames(disease_transform)[1]<-"Tumor_Name"

affinity_cluster_df$SubsetCluster_IDs <- NA

# First round of sub clustering 
for (iter in 1:length(Clusters_Names)) {
  
  print(iter)
  is_large_cluster <- unique(affinity_cluster_df$High_Primary_Cluster_Membership[which(affinity_cluster_df$Cluster_ID==
                                                                                         Clusters_Names[iter])])
  
  if(is_large_cluster=="Yes"){
    
    subset_embedding_df <- as.data.frame(affinity_cluster_df$Tumor_Names[affinity_cluster_df$Cluster_ID==Clusters_Names[iter]])
    colnames(subset_embedding_df)<-"Tumor_Name"
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<- subset_embedding_df %>% dplyr::left_join(disease_transform,by="Tumor_Name")
    rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
    subset_embedding_df<-subset_embedding_df[,c(-1)]
    
    affinity_subset <- apcluster(negDistMat(r=2), subset_embedding_df)
    cat("affinity propogation optimal number of clusters:", length(affinity_subset@clusters), "\n")
    
    subset_affinity_df<-as.data.frame(affinity_subset@idx)
    subset_affinity_df<-as.data.frame(matrix(nrow=1,ncol=2))
    colnames(subset_affinity_df)<-c("Tumor_Names","SubCluster_ID")
    for (iter_subset in 1: length(affinity_subset@clusters)){
      subset_affinity_df[iter_subset,1] <- paste(names(unlist(affinity_subset@clusters[iter_subset])),collapse = "@")
      subset_affinity_df[iter_subset,2] <- iter_subset
    }
    subset_affinity_df<- subset_affinity_df %>% separate_rows(Tumor_Names, sep = '@')
    largest_cluster_id <- max(subset_affinity_df$SubCluster_ID)
    power_of_ten <- floor(log10(largest_cluster_id)) + 1
    subset_affinity_df$SubCluster_ID<-subset_affinity_df$SubCluster_ID/(10^power_of_ten)+Clusters_Names[iter]
    
    for (iter_affinity_cluser in 1: dim(subset_affinity_df)[1]){
      ind_location <- which (affinity_cluster_df$Tumor_Names==subset_affinity_df$Tumor_Names[iter_affinity_cluser])
      affinity_cluster_df$SubsetCluster_IDs[ind_location]<-subset_affinity_df$SubCluster_ID[iter_affinity_cluser]
    }
    
    
  }
}


# Fill up the subclusters columns with NA values. Fill them with primary cluster IDs 
ind_subcluster_na <- which(is.na(affinity_cluster_df$SubsetCluster_IDs),arr.ind = TRUE)

affinity_cluster_df$SubsetCluster_IDs[ind_subcluster_na] <- affinity_cluster_df$Cluster_ID[ind_subcluster_na]

count_subset_table <- as.data.frame(table(affinity_cluster_df$SubsetCluster_IDs))
colnames(count_subset_table) <- c("SubsetCluster_IDs", "Subset_Cluster_Freq")
count_subset_table$SubsetCluster_IDs <- as.numeric(as.character(sub("," , ".", count_subset_table$SubsetCluster_IDs)))

#count_subset_table$SubCluster_IDs<- as.numeric(as.character(count_subset_table$SubsetCluster_IDs))



affinity_cluster_df <- affinity_cluster_df %>% dplyr::left_join(count_subset_table,by ="SubsetCluster_IDs")
affinity_cluster_annotation <- affinity_cluster_df %>% dplyr::select(Tumor_Names,Pediatric_SubsetCluster_ID,SubsetCluster_IDs)


####### CHECK IF NCIT OR WHO 

# Read NCIT Terms and WHO Terms with embedding
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))

NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

affinity_cluster_annotation$NCIT_Tumor<-"No"
affinity_cluster_annotation$WHO_Tumor<-"No"

for (iter in 1:dim(affinity_cluster_annotation)[1]){
  if(affinity_cluster_annotation$Tumor_Names[iter] %in% tolower(NCIT_embedding_df$Disease)){
    affinity_cluster_annotation$NCIT_Tumor[iter] <- "Yes"
  }else if(affinity_cluster_annotation$Tumor_Names[iter] %in% tolower(WHO_embedding_df$Disease)){
    affinity_cluster_annotation$WHO_Tumor[iter] <- "Yes"
  }
}


# write files 
#save(d.apclus2,file = paste(intermediate_dir,"/d.apclus2.RData",sep=""))
save(affinity_cluster_df,file = paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
save(affinity_cluster_annotation,file = paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))

save.image(file = "script5_affinitycluster.RData")
