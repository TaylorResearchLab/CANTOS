# Affinity cluster reassignment after outlier detection
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


affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_df %>% dplyr::select(nct_id,Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_df %>% dplyr::select(nct_id,Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)


outlier_cluster_ada2 <- affinity_cluster_ADA2_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Cluster_ID)%>%distinct()

for(iter in 1:dim(outlier_cluster_ada2)[1]){
  current_cluster_id <- outlier_cluster_ada2$Cluster_ID[iter]
  
  cluster_subset <- affinity_cluster_ADA2_df%>%filter(Cluster_ID==current_cluster_id)
  
  ind_outliers<- which(cluster_subset$Isolation_Outlier=="Yes" | cluster_subset$LOF_Outlier=="Yes")
  
  new_suffix_cluster_id <-1:length(ind_outliers)
  new_cluster_id<- paste(current_cluster_id,new_suffix_cluster_id,sep=";")
  cluster_subset$Cluster_ID[ind_outliers]<-new_cluster_id
  
  cluster_subset<-cluster_subset%>%dplyr::select(Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)
  
  for(iter_subset in ind_outliers){
  cluster_subset$who_cluster_label[iter_subset]<-cluster_subset$WHO_Matches[iter_subset]
  cluster_subset$ncit_cluster_label[iter_subset]<-cluster_subset$NCIT_Matches[iter_subset]
  }
  
  affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df %>% rows_update(cluster_subset,by="Tumor_Names")
  
 
}




cluster_labels<-unique(affinity_cluster_ADA2_reassigned_df$Cluster_ID)
cluster_labels<- as.data.frame(cbind(cluster_labels,c(1:length(cluster_labels))))
colnames(cluster_labels)<-c("Cluster_ID","updated_ID")
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df %>% dplyr::left_join(cluster_labels,by="Cluster_ID")
affinity_cluster_ADA2_reassigned_df<-affinity_cluster_ADA2_reassigned_df[,c(1,2,10,4:9,3)]

affinity_cluster_ADA2_reassigned_df$updated_ID<-as.numeric(affinity_cluster_ADA2_reassigned_df$updated_ID)


#########

outlier_cluster_v3 <- affinity_cluster_v3_df%>%filter(Isolation_Outlier=="Yes" | LOF_Outlier=="Yes")%>%dplyr::select(Cluster_ID)%>%distinct()

for(iter in 1:dim(outlier_cluster_v3)[1]){
  current_cluster_id <- outlier_cluster_v3$Cluster_ID[iter]
  
  cluster_subset <- affinity_cluster_v3_df%>%filter(Cluster_ID==current_cluster_id)
  
  ind_outliers<- which(cluster_subset$Isolation_Outlier=="Yes" | cluster_subset$LOF_Outlier=="Yes")
  
  new_suffix_cluster_id <-1:length(ind_outliers)
  new_cluster_id<- paste(current_cluster_id,new_suffix_cluster_id,sep=";")
  cluster_subset$Cluster_ID[ind_outliers]<-new_cluster_id
  
  cluster_subset<-cluster_subset%>%dplyr::select(Tumor_Names,Cluster_ID,WHO_Matches,WHO_distance,NCIT_Matches,NCIT_distance,who_cluster_label,ncit_cluster_label)
  
  for(iter_subset in ind_outliers){
    cluster_subset$who_cluster_label[iter_subset]<-cluster_subset$WHO_Matches[iter_subset]
    cluster_subset$ncit_cluster_label[iter_subset]<-cluster_subset$NCIT_Matches[iter_subset]
  }

  affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df %>% rows_update(cluster_subset,by="Tumor_Names")
}

cluster_labels<-unique(affinity_cluster_v3_reassigned_df$Cluster_ID)
cluster_labels<- as.data.frame(cbind(cluster_labels,c(1:length(cluster_labels))))
colnames(cluster_labels)<-c("Cluster_ID","updated_ID")
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df %>% dplyr::left_join(cluster_labels,by="Cluster_ID")
affinity_cluster_v3_reassigned_df<-affinity_cluster_v3_reassigned_df[,c(1,2,10,4:9,3)]

affinity_cluster_v3_reassigned_df$updated_ID<-as.numeric(affinity_cluster_v3_reassigned_df$updated_ID)

#####
write.csv(affinity_cluster_ADA2_reassigned_df,paste(intermediate_dir,"/affinity_cluster_ADA2_reassigned_df.csv",sep=""))
write.csv(affinity_cluster_v3_reassigned_df,paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))

