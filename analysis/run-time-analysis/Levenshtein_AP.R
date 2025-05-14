method_name <- "Levenshtein+AP-WHO-ALL"  
start_time <- Sys.time()
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
results_all_dir <- file.path(analysis_dir,"results")
results_5th_dir <- file.path(analysis_dir,"results_5th")


source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))
source(paste(util_dir,"/edit_distance_nested_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))

load(paste(data_dir,"/dissimilarity_matrix_lv.RData",sep=""))
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")
NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])
ct_disease_df <- read_csv(paste(intermediate_dir,"/ct_disease_df.csv",sep=""))
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]
#tumor_nct_map_5th<-read.csv(paste(data_dir,"/tumor_nct_map_5thed.csv",sep=""))
#tumor_nct_map_all<-read.csv(paste(data_dir,"/tumor_nct_map_all.csv",sep=""))


tumor_names_5th<- unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_5th$Tumor_Names))




dissimilarity_matrix_lv_5th<- dissimilarity_matrix_lv %>% dplyr::filter(rownames(dissimilarity_matrix_lv) %in% 
                                                                          tumor_names_5th )
dissimilarity_matrix_lv_5th<- dissimilarity_matrix_lv_5th %>%dplyr::select(any_of(tumor_names_5th))
dissimilarity_matrix_lv_5th<-as.matrix(dissimilarity_matrix_lv_5th)


dissimilarity_matrix_lv<-as.matrix(dissimilarity_matrix_lv)


simmilarity_matrix_lv <- 1- dissimilarity_matrix_lv
simmilarity_matrix_lv_5th <- 1- dissimilarity_matrix_lv_5th 
#######################

print("computing with WHO all")
apclust_lv_all <- apcluster(simmilarity_matrix_lv) #11:19 am -2:44 pm 
affinity_cluster_lv_all_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_lv_all_df)<-c("Tumor_Names","Cluster_ID")

start_time <- Sys.time()
for (iter in 1: length(apclust_lv_all@clusters)){
  affinity_cluster_lv_all_df[iter,1] <- paste(names(unlist(apclust_lv_all@clusters[iter])),collapse = "@")
  affinity_cluster_lv_all_df[iter,2] <- iter
}
affinity_cluster_lv_all_df<- affinity_cluster_lv_all_df %>% separate_rows(Tumor_Names, sep = '@')


nested_affinity_cluster_lv<- edit_distance_nested_cluster(affinity_cluster_lv_all_df,simmilarity_matrix_lv)
nested_affinity_cluster_lv<-compute_silhouette(cluster_df = nested_affinity_cluster_lv,dist_mat = dissimilarity_matrix_lv)
mean_freq_lv <- nested_affinity_cluster_lv %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(mean_freq_lv,by="Cluster_ID")
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(7,1:6)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(-7)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv %>% dplyr::select(Tumor_Names,Cluster_ID)

dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
dissimilarity_matrix_lv_who_all <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Terms_All$Tumor_Names))
dissimilarity_matrix_lv_ncit <- dissimilarity_matrix_lv %>% dplyr:::select(one_of(NCIT_Terms))

index_min_who_all_lv <- as.matrix(apply(dissimilarity_matrix_lv_who_all, 1, which.min))
who_all_match_lv_df <- cbind(rownames(dissimilarity_matrix_lv_who_all))
colnames(who_all_match_lv_df)<-"Tumor_Names"
who_all_match_lv_df <-as.data.frame(who_all_match_lv_df)
who_all_match_lv_df$WHO_Matches<- NA
who_all_match_lv_df$WHO_distance<-NA

for (iter in 1: dim(who_all_match_lv_df)[1]){
  
  who_all_match_lv_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_lv_who_all)[index_min_who_all_lv[iter]]
  who_all_match_lv_df$WHO_distance[iter]<-dissimilarity_matrix_lv_who_all[iter,index_min_who_all_lv[iter]]
  
}

index_min_ncit_lv <- as.matrix(apply(dissimilarity_matrix_lv_ncit, 1, which.min))
ncit_match_lv_df <- cbind(rownames(dissimilarity_matrix_lv_ncit))
colnames(ncit_match_lv_df)<-"Tumor_Names"
ncit_match_lv_df <-as.data.frame(ncit_match_lv_df)
ncit_match_lv_df$NCIT_Matches<- NA
ncit_match_lv_df$NCIT_distance<-NA

for (iter in 1: dim(ncit_match_lv_df)[1]){
  
  ncit_match_lv_df$NCIT_Matches[iter] <- colnames(dissimilarity_matrix_lv_ncit)[index_min_ncit_lv[iter]]
  ncit_match_lv_df$NCIT_distance[iter]<-dissimilarity_matrix_lv_ncit[iter,index_min_ncit_lv[iter]]
  
}

nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(who_all_match_lv_df,by="Tumor_Names")
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(ncit_match_lv_df,by="Tumor_Names")

nested_affinity_cluster_lv<- cluster_label_assignment_refined(nested_affinity_cluster_lv)
nested_affinity_cluster_lv<-outlier_detection_edit_dist(nested_affinity_cluster_lv,dissimilarity_matrix_lv)
nested_affinity_cluster_lv_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_lv)
nested_affinity_cluster_lv_reassigned_short <- nested_affinity_cluster_lv_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label,ncit_cluster_label)


end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
print(elapsed_time)


benchmark <- data.frame(
  Method = method_name,
  Runtime_sec = round(elapsed_time, 2),
  stringsAsFactors = FALSE
)

output_file <- paste(analysis_dir,"/run-time-analysis/runtime_benchmarks_all_methods.csv",sep="")
if (!file.exists(output_file)) {
  write.csv(benchmark, output_file, row.names = FALSE)
} else {
  write.table(benchmark, output_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}




