method_name <- "cosine+AP-WHO-ALL"  
start_time_1 <- Sys.time()
# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(apcluster)
  library(stringdist)
  library(tidyverse)
  library(magrittr)
  library(isotree)
  library(dbscan)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")


source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))
source(paste(util_dir,"/edit_distance_nested_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))

load(paste(data_dir,"/dissimilarity_matrix_cosine.RData",sep=""))
dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")
NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])
ct_disease_df <- read_csv(paste(intermediate_dir,"/ct_disease_df.csv",sep=""))
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]


dissimilarity_matrix_cosine<-as.matrix(dissimilarity_matrix_cosine)


simmilarity_matrix_cosine <- 1- dissimilarity_matrix_cosine
#######################

print("computing with WHO all")
apclust_cosine_all <- apcluster(simmilarity_matrix_cosine) #11:19 am -2:44 pm 
affinity_cluster_cosine_all_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_cosine_all_df)<-c("Tumor_Names","Cluster_ID")

for (iter in 1: length(apclust_cosine_all@clusters)){
  affinity_cluster_cosine_all_df[iter,1] <- paste(names(unlist(apclust_cosine_all@clusters[iter])),collapse = "@")
  affinity_cluster_cosine_all_df[iter,2] <- iter
}
affinity_cluster_cosine_all_df<- affinity_cluster_cosine_all_df %>% separate_rows(Tumor_Names, sep = '@')

end_time_1 <- Sys.time()
elapsed_time_1 <- as.numeric(difftime(end_time_1, start_time_1, units = "secs"))

save.image(paste(run_time_analysis,"/cosine_AP/cosine_AP_ALL.RData",sep=""))


start_time_2 <- Sys.time()


nested_affinity_cluster_cosine<- edit_distance_nested_cluster(affinity_cluster_cosine_all_df,simmilarity_matrix_cosine)
nested_affinity_cluster_cosine<-compute_silhouette(cluster_df = nested_affinity_cluster_cosine,dist_mat = dissimilarity_matrix_cosine)
mean_freq_cosine <- nested_affinity_cluster_cosine %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(mean_freq_cosine,by="Cluster_ID")
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(7,1:6)]
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(-7)]
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine %>% dplyr::select(Tumor_Names,Cluster_ID)


end_time_2 <- Sys.time()
save.image(paste(run_time_analysis,"/cosine_AP/cosine_AP_ALL.RData",sep=""))

elapsed_time_2 <- as.numeric(difftime(end_time_2, start_time_2, units = "secs"))

start_time3<-Sys.time()

dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_cosine_who_all <- dissimilarity_matrix_cosine %>% dplyr::select(one_of(WHO_Terms_All$Tumor_Names))
dissimilarity_matrix_cosine_ncit <- dissimilarity_matrix_cosine %>% dplyr:::select(one_of(NCIT_Terms))

index_min_who_all_cosine <- as.matrix(apply(dissimilarity_matrix_cosine_who_all, 1, which.min))
who_all_match_cosine_df <- cbind(rownames(dissimilarity_matrix_cosine_who_all))
colnames(who_all_match_cosine_df)<-"Tumor_Names"
who_all_match_cosine_df <-as.data.frame(who_all_match_cosine_df)
who_all_match_cosine_df$WHO_Matches<- NA
who_all_match_cosine_df$WHO_distance<-NA

for (iter in 1: dim(who_all_match_cosine_df)[1]){
  
  who_all_match_cosine_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_cosine_who_all)[index_min_who_all_cosine[iter]]
  who_all_match_cosine_df$WHO_distance[iter]<-dissimilarity_matrix_cosine_who_all[iter,index_min_who_all_cosine[iter]]
  
}

index_min_ncit_cosine <- as.matrix(apply(dissimilarity_matrix_cosine_ncit, 1, which.min))
ncit_match_cosine_df <- cbind(rownames(dissimilarity_matrix_cosine_ncit))
colnames(ncit_match_cosine_df)<-"Tumor_Names"
ncit_match_cosine_df <-as.data.frame(ncit_match_cosine_df)
ncit_match_cosine_df$NCIT_Matches<- NA
ncit_match_cosine_df$NCIT_distance<-NA

for (iter in 1: dim(ncit_match_cosine_df)[1]){
  
  ncit_match_cosine_df$NCIT_Matches[iter] <- colnames(dissimilarity_matrix_cosine_ncit)[index_min_ncit_cosine[iter]]
  ncit_match_cosine_df$NCIT_distance[iter]<-dissimilarity_matrix_cosine_ncit[iter,index_min_ncit_cosine[iter]]
  
}

nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(who_all_match_cosine_df,by="Tumor_Names")
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(ncit_match_cosine_df,by="Tumor_Names")

nested_affinity_cluster_cosine<- cluster_label_assignment_refined(nested_affinity_cluster_cosine)
nested_affinity_cluster_cosine<-outlier_detection_edit_dist(nested_affinity_cluster_cosine,dissimilarity_matrix_cosine)
nested_affinity_cluster_cosine_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_cosine)
nested_affinity_cluster_cosine_reassigned_short <- nested_affinity_cluster_cosine_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label,ncit_cluster_label)


end_time_3 <- Sys.time()
save.image(paste(run_time_analysis,"/cosine_AP/cosine_AP_ALL.RData",sep=""))

elapsed_time_3 <- as.numeric(difftime(end_time_3, start_time3, units = "secs"))
elapsed_time<-elapsed_time_1+elapsed_time_2+elapsed_time_3

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

Runtime_cosine_AP_all = round(elapsed_time, 2)
Memory_gb_cosine_AP_all=round(sum(sizes) / (1024^2), 3)/1000

save.image(paste(run_time_analysis,"/cosine_AP/cosine_AP_ALL.RData",sep=""))
