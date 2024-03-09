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
results_dir <- file.path(analysis_dir,"results")
# Load affinity Cluster
#load(paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))
source(paste(util_dir,"/nested_affinity_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment.R",sep=""))



########################################*************************************###########
# Load PCA Embeddings of CT , WHO, NCIT
disease_transform<- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )
colnames(disease_transform)[1]<-"Tumor_Name"
rownames(disease_transform)<-disease_transform$Tumor_Name # Needed for AP Clust


# Set Seed
set.seed(13)
#affinity cluster
d.apclus2 <- apcluster(negDistMat(r=2), disease_transform) # 1 hr 28 mins 11:28 pm - 12:08 pm

dist_euclidean<- dist(disease_transform,method = "euclidean")
dist_euclidean<-as.matrix(dist_euclidean)
simmilarity_euclidean<- 1/(1+dist_euclidean)
af_clust_euclidean <- apcluster(simmilarity_euclidean) # 1:24 am - 2:55 am still going....
cat("affinity propogation optimal number of clusters:", length(af_clust_euclidean@clusters), "\n")

cat("affinity propogation optimal number of clusters:", length(d.apclus2@clusters), "\n") #1113 clusters 

affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(d.apclus2@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(d.apclus2@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
affinity_cluster_df$Cluster_ID<-as.character(affinity_cluster_df$Cluster_ID)


affinity_cluster_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(af_clust_euclidean@clusters)){
  affinity_cluster_df[iter,1] <- paste(names(unlist(af_clust_euclidean@clusters[iter])),collapse = "@")
  affinity_cluster_df[iter,2] <- iter
}
affinity_cluster_df<- affinity_cluster_df %>% separate_rows(Tumor_Names, sep = '@')
affinity_cluster_df$Cluster_ID<-as.character(affinity_cluster_df$Cluster_ID)



# Find cluster membership frequencies 
#affinity_cluster_df$Cluster_Total_Members <- NA
cluster_frequency_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
z_scores<- (cluster_frequency_table$Primary_Cluster_Frequency-mean(cluster_frequency_table$Primary_Cluster_Frequency))/sd(cluster_frequency_table$Primary_Cluster_Frequency)
cluster_frequency_table$z_scores<-z_scores

ind_min_zscore<- which(cluster_frequency_table$z_scores < 2.5)
max_cluster_member <- max(cluster_frequency_table$Primary_Cluster_Frequency[ind_min_zscore])

#median_cluster_frequency <- median(cluster_frequency_table$Primary_Cluster_Frequency)

#large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>median_cluster_frequency)]
large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]

converge_list<-list()

while(length(large_cluster_labels)>0){
  print(length(large_cluster_labels))
for(iter in 1:length(large_cluster_labels)){
  Clusters_Names=large_cluster_labels[iter]
  subset_embedding_df <- as.data.frame(affinity_cluster_df$Tumor_Names[affinity_cluster_df$Cluster_ID==Clusters_Names])
  colnames(subset_embedding_df)<-"Tumor_Name"
  rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
  subset_embedding_df<- subset_embedding_df %>% dplyr::left_join(disease_transform,by="Tumor_Name")
  rownames(subset_embedding_df)<-subset_embedding_df$Tumor_Name
  subset_embedding_df<-subset_embedding_df[,c(-1)]
  
  result_run_aff<-run_affinity_clustering(Clusters_Names,subset_embedding_df)
  
  flag_converge <- result_run_aff[[1]]
  subset_affinity_df<-result_run_aff[[2]]
  
  if(flag_converge=="No"){
  for (iter_nested_affinity_cluser in 1: dim(subset_affinity_df)[1]){
      ind_location <- which (affinity_cluster_df$Tumor_Names==subset_affinity_df$Tumor_Names[iter_nested_affinity_cluser])
      affinity_cluster_df$Cluster_ID[ind_location]<-subset_affinity_df$SubCluster_ID[iter_nested_affinity_cluser]
   }
  }else if(flag_converge=="Yes"){
    converge_list<-append(Clusters_Names,converge_list)
  }


 }
  cluster_frequency_table <- as.data.frame(table(affinity_cluster_df$Cluster_ID))
  colnames(cluster_frequency_table)<- c("Cluster_ID","Primary_Cluster_Frequency")
  cluster_frequency_table$Cluster_ID<-as.character(cluster_frequency_table$Cluster_ID)
  large_cluster_labels<- cluster_frequency_table$Cluster_ID[which(cluster_frequency_table$Primary_Cluster_Frequency>max_cluster_member)]
  large_cluster_labels<-setdiff(large_cluster_labels, unlist(converge_list))
}
  
lymphoma_leukemia_strings <- c("leukemia", "lymphoma", "leukemias", "lymphomas", "leukaemia", "leukaemias",
                               "leuk")

Lymphoma_Lukemia_Clusters <- affinity_cluster_df %>% dplyr::filter(str_detect(Tumor_Names,paste(strings, collapse = "|")))

Lymphoma_Lukemia_Clusters_Labels <- unique(Lymphoma_Lukemia_Clusters$Cluster_ID)

Lymphoma_Lukemia_Clusters <- affinity_cluster_df %>% dplyr::filter(Cluster_ID %in% Lymphoma_Lukemia_Clusters_Labels )
  
  
  # largest_cluster_id <- max(subset_affinity_df$SubCluster_ID)
  # power_of_ten <- floor(log10(largest_cluster_id)) + 1
  # subset_affinity_df$SubCluster_ID<-subset_affinity_df$SubCluster_ID/(10^power_of_ten)+Clusters_Names[iter]
  
  for (iter_affinity_cluser in 1: dim(subset_affinity_df)[1]){
    ind_location <- which (affinity_cluster_df$Tumor_Names==subset_affinity_df$Tumor_Names[iter_affinity_cluser])
    affinity_cluster_df$SubsetCluster_IDs[ind_location]<-subset_affinity_df$SubCluster_ID[iter_affinity_cluser]
  }
  
  
  



# # Find if tumor has children and pediatric terms in their tumor names and make a seperate column for pediatric cluster
# ind_ped <- c(which(str_detect(affinity_cluster_df$Tumor_Names, "childhood")),
#              which(str_detect(affinity_cluster_df$Tumor_Names, "children")),
#              which(str_detect(affinity_cluster_df$Tumor_Names, "child")),
#              which(str_detect(affinity_cluster_df$Tumor_Names, "pediatric")),
#              which(str_detect(affinity_cluster_df$Tumor_Names, "paediatric")))
# 
# affinity_cluster_df$contains_pediatric_string <- "No"
# 
# affinity_cluster_df$contains_pediatric_string[ind_ped]<-"Yes"
# 
# pediatric_cluster_ID<-0
# 
# affinity_cluster_df$Pediatric_SubsetCluster_ID <- NA  
# affinity_cluster_df$Pediatric_SubsetCluster_ID[ind_ped]<-pediatric_cluster_ID


#Nested affinity clusters
# barplot(height=affinity_cluster_df$Primary_Cluster_Frequency, names=affinity_cluster_df$Cluster_ID, 
#         xlab="cluster_id", 
#         ylab="frequency", 
#         main="Cluster_Frequency", 
#         ylim=c(0,50)
# )

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











########################################*************************************###########
# load disease_transfor
disease_transform <- read.csv(paste(intermediate_dir,"/disease_transform_pca.csv",sep=""))
colnames(disease_transform)[1]<-"Tumor_Name"

# Load data
ncit_match_df <- read.csv(paste(intermediate_dir,"/ncit_match_df.csv",sep=""))
who_match_df <- read.csv(paste(intermediate_dir,"/who_ct_distance_mat.csv",sep=""))


# Need for second run of clustering

affinity_cluster_nested <- affinity_cluster_annotation %>% dplyr::select(Tumor_Names,Pediatric_SubsetCluster_ID, SubsetCluster_IDs)
affinity_cluster_nested <- nested_affinity_cluster(n=3,affinity_cluster_nested,disease_transform)


# Silos scores Affinity Cluster


affinity_cluster_nested<- affinity_cluster_nested %>% dplyr::left_join(ncit_match_df,by="Tumor_Names")
affinity_cluster_nested <- affinity_cluster_nested %>%dplyr::left_join(who_match_df,by="Tumor_Names")



affinity_cluster_nested <- affinity_cluster_nested %>% dplyr::mutate(assigned_class = case_when(ncit_distance < WHO_distance ~ NCIT_Matches,
                                                                                                ncit_distance > WHO_distance ~ WHO_Matches,
                                                                                                TRUE ~ "Both"))

# Cluster voting
affinity_cluster_nested<- cluster_label_assignment(affinity_cluster_nested)

disease_affinity_cluster_table<- affinity_cluster_nested %>% dplyr::select(Tumor_Names,cluster_label)

# Write
write.csv(affinity_cluster_nested,paste(intermediate_dir,"/affinity_cluster_nested.csv",sep=""))

save.image(file = "script6_affinitycluster.RData")
# write files 
#save(d.apclus2,file = paste(intermediate_dir,"/d.apclus2.RData",sep=""))
save(affinity_cluster_df,file = paste(intermediate_dir,"/affinity_cluster_df.RData",sep=""))
save(affinity_cluster_annotation,file = paste(intermediate_dir,"/affinity_cluster_annotation.RData",sep=""))

# Silos not computed
source("~/Desktop/MTP_Paper/CT-Embedding-Paper/util/compute_silhouette.R")
affinity_cluster_df2<-affinity_cluster_df
colnames(affinity_cluster_df2)[2]<-"SubsetCluster_IDs"
affinity_cluster_df2<-compute_silhouette(affinity_cluster_df2,dist_euclidean) # Change colname to sublu
save.image(file = "script6_affinitycluster.RData")
mean_freq_af <- affinity_cluster_df2 %>%dplyr::group_by(SubsetCluster_IDs) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
affinity_cluster_df2<- affinity_cluster_df2 %>% dplyr::left_join(mean_freq_af,by="SubsetCluster_IDs")

benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

cluster_ind_benchmark_tumor <- affinity_cluster_df2$SubsetCluster_IDs[affinity_cluster_df2$Tumor_Names %in% benchmark_tumors]

display_table_benchmark_af <- affinity_cluster_df2 %>% filter(SubsetCluster_IDs %in% cluster_ind_benchmark_tumor)
display_table_benchmark_af<- display_table_benchmark_af[order(display_table_benchmark_af$SubsetCluster_IDs),]
rownames(display_table_benchmark_af)<-NULL

write.csv(display_table_benchmark_af,paste(results_dir,"/display_table_benchmark_embedding_af.csv",sep=""))