suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
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

# Load data 
affinity_cluster_ADA2_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_ADA2_df.csv",sep=""))
affinity_cluster_v3_df <- read.csv(paste(intermediate_dir,"/affinity_cluster_v3_df.csv",sep=""))

affinity_cluster_ADA2_df<-affinity_cluster_ADA2_df[,c(-1)]
affinity_cluster_v3_df<- affinity_cluster_v3_df[,c(-1)]


affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_df %>% dplyr::select(Tumor_Names,Cluster_ID,assigned_class,suggested_cluster_label)
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_df %>% dplyr::select(Tumor_Names,Cluster_ID,assigned_class,suggested_cluster_label)

# Find all the tumors with outlier
# outlier_tumor_ada2 <- affinity_cluster_ADA2_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Tumor_Names)
# outlier_tumor_v3 <- affinity_cluster_v3_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Tumor_Names)


outlier_cluster_ada2 <- affinity_cluster_ADA2_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Cluster_ID)%>%distinct()
outlier_cluster_v3 <- affinity_cluster_v3_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Cluster_ID)%>%distinct()




for(iter in 1:dim(outlier_cluster_ada2)[1]){
  current_cluster_id <- outlier_cluster_ada2$Cluster_ID[iter]
  
  cluster_subset <- affinity_cluster_ADA2_df%>%filter(Cluster_ID==current_cluster_id)
  
  ind_outliers<- which(cluster_subset$Isolation_Outlier=="Yes" | cluster_subset$LOF_Outlier=="Yes")
  
  new_suffix_cluster_id <-1:length(ind_outliers)
  new_cluster_id<- paste(current_cluster_id,new_suffix_cluster_id,sep=";")
  cluster_subset$Cluster_ID[ind_outliers]<-new_cluster_id
  
  cluster_subset<-cluster_subset%>%dplyr::select(Tumor_Names,Cluster_ID,assigned_class,suggested_cluster_label)
  
  for(iter_subset in 1:dim(cluster_subset)[1]){
      cluster_subset$suggested_cluster_label[iter_subset]<-cluster_subset$assigned_class[iter_subset]
  }
  
  table_frequency_assigned_class<- as.data.frame(table(cluster_subset$assigned_class))
  table_frequency_assigned_class$Var1<-as.character(table_frequency_assigned_class$Var1)
  high_frequency_assigned_class<- table_frequency_assigned_class$Var1[which(table_frequency_assigned_class$Freq>1)]
  
  if(length(high_frequency_assigned_class)>0){
    for(iter_high_frequency in 1:length(high_frequency_assigned_class)){
      ind_matches_sub <- which(cluster_subset$assigned_class==high_frequency_assigned_class[iter_high_frequency]) 
      cluster_subset$Cluster_ID[ind_matches_sub]<-cluster_subset$Cluster_ID[ind_matches_sub[1]]
      cluster_subset$suggested_cluster_label[ind_matches_sub]<-cluster_subset$suggested_cluster_label[ind_matches_sub[1]]
    }
  }
  affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df %>% rows_update(cluster_subset,by="Tumor_Names")
}

cluster_labels<-unique(affinity_cluster_ADA2_reassigned_df$Cluster_ID)
cluster_labels<- as.data.frame(cbind(cluster_labels,c(1:length(cluster_labels))))
colnames(cluster_labels)<-c("Cluster_ID","updated_ID")
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df %>% dplyr::left_join(cluster_labels,by="Cluster_ID")
for(iter in 1:length(cluster_labels)){
  ind_location_label_ADA2 <- which(affinity_cluster_ADA2_reassigned_df$Cluster_ID==cluster_labels[iter])
  affinity_cluster_ADA2_reassigned_df$Cluster_ID[ind_location_label_ADA2]<-iter
}






