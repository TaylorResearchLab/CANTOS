method_name <- "Levenshtein+AP-WHO-5th"  
# Load libraries
start_time <- Sys.time()
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
results_dir <- file.path(analysis_dir,"results_5th")
plots_dir <- file.path(root_dir,"plots")
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")


source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))
source(paste(util_dir,"/edit_distance_nested_cluster.R",sep=""))
source(paste(util_dir,"/cluster_label_assignment_refined.R",sep=""))
source(paste(util_dir,"/outlier_detection_edit_dist.R",sep=""))
source(paste(util_dir,"/edit_distance_cluster_reassignment.R",sep=""))

print("LV_5th")

# Load data
load(paste(data_dir,"/dissimilarity_matrix_lv.RData",sep=""))
#load(paste(data_dir,"/dissimilarity_matrix_cosine.RData",sep=""))
#load(paste(data_dir,"/dissimilarity_matrix_jw.RData",sep=""))
# load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_lv.RData")
# load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_cosine.RData")
# load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_jw.RData")
# 

# Convert dissimilarity matrix to only contain who 5th edition
#dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
#dissimilarity_matrix_jw<-as.data.frame(dissimilarity_matrix_jw)


WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])
ct_disease_df <- read_csv(paste(intermediate_dir,"/ct_disease_df.csv",sep=""))
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]

tumor_names<- unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_5th$Tumor_Names))

dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>% dplyr::filter(rownames(dissimilarity_matrix_lv) %in% 
                                                                      tumor_names )



dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>%dplyr::select(any_of(tumor_names))

dissimilarity_matrix_lv<-as.matrix(dissimilarity_matrix_lv)



# Compute Similarity matrix for each edit distance
simmilarity_matrix_lv <- 1- dissimilarity_matrix_lv

print("LV_5th APclust")

######### Cluster with LV ########
apclust_lv <- apcluster(simmilarity_matrix_lv) 
affinity_cluster_lv_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_lv_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_lv@clusters)){
  affinity_cluster_lv_df[iter,1] <- paste(names(unlist(apclust_lv@clusters[iter])),collapse = "@")
  affinity_cluster_lv_df[iter,2] <- iter
}
affinity_cluster_lv_df<- affinity_cluster_lv_df %>% separate_rows(Tumor_Names, sep = '@')



end_time <- Sys.time()

elapsed_time_1 <- as.numeric(difftime(end_time, start_time, units = "secs"))


save.image(paste(run_time_analysis,"/LV_AP/LV_AP_5th.RData",sep=""))

start_time2 <- Sys.time()

nested_affinity_cluster_lv<- edit_distance_nested_cluster(affinity_cluster_lv_df,simmilarity_matrix_lv)
nested_affinity_cluster_lv<-compute_silhouette(cluster_df = nested_affinity_cluster_lv,dist_mat = dissimilarity_matrix_lv)
mean_freq_lv <- nested_affinity_cluster_lv %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(mean_freq_lv,by="Cluster_ID")
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(7,1:6)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(-7)]
nested_affinity_cluster_lv<-nested_affinity_cluster_lv %>% dplyr::select(Tumor_Names,Cluster_ID)
end_time2 <- Sys.time()

elapsed_time_2 <- as.numeric(difftime(end_time2, start_time2, units = "secs"))

save.image(paste(run_time_analysis,"/LV_AP/LV_AP_5th.RData",sep=""))


start_time3<-Sys.time()
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
dissimilarity_matrix_lv_who_5th <- dissimilarity_matrix_lv %>% dplyr::select(one_of(WHO_Terms_5th$Tumor_Names))
dissimilarity_matrix_lv_ncit <- dissimilarity_matrix_lv %>% dplyr:::select(one_of(NCIT_Terms))

index_min_who_5th_lv <- as.matrix(apply(dissimilarity_matrix_lv_who_5th, 1, which.min))
who_5th_match_lv_df <- cbind(rownames(dissimilarity_matrix_lv_who_5th))
colnames(who_5th_match_lv_df)<-"Tumor_Names"
who_5th_match_lv_df <-as.data.frame(who_5th_match_lv_df)
who_5th_match_lv_df$WHO_Matches<- NA
who_5th_match_lv_df$WHO_distance<-NA

for (iter in 1: dim(who_5th_match_lv_df)[1]){
  
  who_5th_match_lv_df$WHO_Matches[iter] <- colnames(dissimilarity_matrix_lv_who_5th)[index_min_who_5th_lv[iter]]
  who_5th_match_lv_df$WHO_distance[iter]<-dissimilarity_matrix_lv_who_5th[iter,index_min_who_5th_lv[iter]]
  
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

end_time3 <- Sys.time()
elapsed_time_3 <- as.numeric(difftime(end_time3, start_time3, units = "secs"))


save.image(paste(run_time_analysis,"/LV_AP/LV_AP_5th.RData",sep=""))


start_time4<-Sys.time()



nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(who_5th_match_lv_df,by="Tumor_Names")
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(ncit_match_lv_df,by="Tumor_Names")

nested_affinity_cluster_lv<- cluster_label_assignment_refined(nested_affinity_cluster_lv)
nested_affinity_cluster_lv<-outlier_detection_edit_dist(nested_affinity_cluster_lv,dissimilarity_matrix_lv)
nested_affinity_cluster_lv_reassigned<-edit_distance_cluster_reassignment(nested_affinity_cluster_lv)
nested_affinity_cluster_lv_reassigned_short <- nested_affinity_cluster_lv_reassigned %>% dplyr::select(Tumor_Names,who_cluster_label,ncit_cluster_label)


end_time4 <- Sys.time()


elapsed_time_4 <- as.numeric(difftime(end_time4, start_time4, units = "secs"))

elapsed_time<-elapsed_time_1+elapsed_time_2+elapsed_time_3+elapsed_time_4


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

Runtime_LV_AP_5th = round(elapsed_time, 2)
Memory_gb_LV_AP_5th=round(sum(sizes) / (1024^2), 3)/1000



save.image(paste(run_time_analysis,"/LV_AP/LV_AP_5th.RData",sep=""))









