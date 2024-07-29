# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(doParallel)
  library(foreach)
  library(apcluster)
  library(stringdist)
  library(tidyverse)
  library(magrittr)
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
source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))
source(paste(util_dir,"/edit_distance_nested_cluster.R",sep=""))



# Load data

load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_lv.RData")
load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_cosine.RData")
load("/Users/lahiria/Desktop/MTP_Paper/temp/CT-Large-File-June21/dissimilarity_matrix_jw.RData")


# Convert dissimilarity matrix to only contain who 5th edition
dissimilarity_matrix_cosine<-as.data.frame(dissimilarity_matrix_cosine)
dissimilarity_matrix_lv<-as.data.frame(dissimilarity_matrix_lv)
dissimilarity_matrix_jw<-as.data.frame(dissimilarity_matrix_jw)


WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

NCIT_Terms <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))[,1]
NCIT_Terms<-tolower(NCIT_Terms[c(-1)])

ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_adult_ped.csv",sep=""))
ct_tumor_df<-ct_disease_annot_adult_ped_df%>%filter(validated_cancer_tumor=="Yes")
ct_tumor_df<-ct_tumor_df[,c(-1)]

tumor_names_all<- unique(c(ct_tumor_df$diseases,NCIT_Terms,WHO_Terms_5th$Tumor_Names))

dissimilarity_matrix_cosine<- dissimilarity_matrix_cosine %>% dplyr::filter(rownames(dissimilarity_matrix_cosine) %in% 
                                                                              tumor_names_all )
dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>% dplyr::filter(rownames(dissimilarity_matrix_lv) %in% 
                                                                      tumor_names_all )
dissimilarity_matrix_jw<- dissimilarity_matrix_jw %>% dplyr::filter(rownames(dissimilarity_matrix_jw) %in% 
                                                                      tumor_names_all )


dissimilarity_matrix_cosine<- dissimilarity_matrix_cosine %>%dplyr::select(any_of(tumor_names_all))
dissimilarity_matrix_lv<- dissimilarity_matrix_lv %>%dplyr::select(any_of(tumor_names_all))
dissimilarity_matrix_jw<- dissimilarity_matrix_jw %>% dplyr::select(any_of(tumor_names_all))

dissimilarity_matrix_cosine<-as.matrix(dissimilarity_matrix_cosine)
dissimilarity_matrix_lv<-as.matrix(dissimilarity_matrix_lv)
dissimilarity_matrix_jw<-as.matrix(dissimilarity_matrix_jw)

# 
cluster_results_lv<- read.csv(paste(results_dir,"/cluster_lv.csv",sep=""))
cluster_results_lv<-cluster_results_lv[,c(-1)]
cluster_results_jw<- read.csv(paste(results_dir,"/cluster_jw.csv",sep=""))
cluster_results_jw<-cluster_results_jw[,c(-1)]
cluster_results_cosine<- read.csv(paste(results_dir,"/cluster_cosine.csv",sep=""))
cluster_results_cosine<-cluster_results_cosine[,c(-1)]

# Compute Similarity matrix for each edit distance
simmilarity_matrix_cosine <- 1 - dissimilarity_matrix_cosine
simmilarity_matrix_jw <- 1-dissimilarity_matrix_jw
simmilarity_matrix_lv <- 1- dissimilarity_matrix_lv




######### Cluster with LV ########
apclust_lv <- apcluster(simmilarity_matrix_lv) #11:19 am -2:44 pm 
affinity_cluster_lv_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_lv_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_lv@clusters)){
  affinity_cluster_lv_df[iter,1] <- paste(names(unlist(apclust_lv@clusters[iter])),collapse = "@")
  affinity_cluster_lv_df[iter,2] <- iter
}
affinity_cluster_lv_df<- affinity_cluster_lv_df %>% separate_rows(Tumor_Names, sep = '@')



######### Cluster with jw ########
apclust_jw <- apcluster(simmilarity_matrix_jw) #2:53 pm 
affinity_cluster_jw_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_jw_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_jw@clusters)){
  affinity_cluster_jw_df[iter,1] <- paste(names(unlist(apclust_jw@clusters[iter])),collapse = "@")
  affinity_cluster_jw_df[iter,2] <- iter
}
affinity_cluster_jw_df<- affinity_cluster_jw_df %>% separate_rows(Tumor_Names, sep = '@')

######### Cluster with cosine ########
apclust_cosine <- apcluster(simmilarity_matrix_cosine)#11:07 am - 1:56 pm
affinity_cluster_cosine_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_cosine_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_cosine@clusters)){
  affinity_cluster_cosine_df[iter,1] <- paste(names(unlist(apclust_cosine@clusters[iter])),collapse = "@")
  affinity_cluster_cosine_df[iter,2] <- iter
}
affinity_cluster_cosine_df<- affinity_cluster_cosine_df %>% separate_rows(Tumor_Names, sep = '@')

# nested clustering based on size of cluster
nested_affinity_cluster_lv<- edit_distance_nested_cluster(affinity_cluster_lv_df,simmilarity_matrix_lv)
nested_affinity_cluster_jw<- edit_distance_nested_cluster(affinity_cluster_jw_df,simmilarity_matrix_jw)
nested_affinity_cluster_cosine<- edit_distance_nested_cluster(affinity_cluster_cosine_df,simmilarity_matrix_cosine)


# Compute Silo Dist

nested_affinity_cluster_lv<-compute_silhouette(cluster_df = nested_affinity_cluster_lv,dist_mat = dissimilarity_matrix_lv)
nested_affinity_cluster_jw<-compute_silhouette(cluster_df = nested_affinity_cluster_jw,dist_mat = dissimilarity_matrix_jw)
nested_affinity_cluster_cosine<-compute_silhouette(cluster_df = nested_affinity_cluster_cosine,dist_mat = dissimilarity_matrix_cosine)


mean_freq_lv <- nested_affinity_cluster_lv %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
mean_freq_jw <- nested_affinity_cluster_jw %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n())
mean_freq_cosine <- nested_affinity_cluster_cosine %>% dplyr::select(Cluster_ID, silhouette_score) %>% dplyr::group_by(Cluster_ID) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n())


nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(mean_freq_lv,by="Cluster_ID")
nested_affinity_cluster_jw<- nested_affinity_cluster_jw %>% dplyr::left_join(mean_freq_jw,by="Cluster_ID")
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(mean_freq_cosine,by="Cluster_ID")


ct_disease_df <- read_csv(paste(intermediate_dir,"/ct_disease_df.csv",sep=""))



#add NCT_IDs
nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))
nested_affinity_cluster_jw<- nested_affinity_cluster_jw %>%  left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>%   left_join(ct_disease_df, by=c("Tumor_Names"="diseases"))

nested_affinity_cluster_lv<-nested_affinity_cluster_lv[,c(7,1:6)]
nested_affinity_cluster_jw<-nested_affinity_cluster_jw[,c(7,1:6)]
nested_affinity_cluster_cosine<-nested_affinity_cluster_cosine[,c(7,1:6)]

# Select benchmarks 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

subcluster_lv <- nested_affinity_cluster_lv$Cluster_ID[nested_affinity_cluster_lv$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_lv <- nested_affinity_cluster_lv %>% dplyr::filter(Cluster_ID %in% subcluster_lv)

subcluster_jw <- nested_affinity_cluster_jw$Cluster_ID[nested_affinity_cluster_jw$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_jw <- nested_affinity_cluster_jw %>% dplyr::filter(Cluster_ID %in% subcluster_jw)

subcluster_cosine <- nested_affinity_cluster_cosine$Cluster_ID[nested_affinity_cluster_cosine$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_cosine <- nested_affinity_cluster_cosine %>% dplyr::filter(Cluster_ID %in% subcluster_cosine)


# Plot 

clust_plot_lv <- benchmark_aff_clust_lv %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) 
clust_plot_lv<-unique(clust_plot_lv)
p_lv <- ggplot(clust_plot_lv, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_lv$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Levenshtein distance clusters")


clust_plot_jw <- benchmark_aff_clust_jw %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) 
clust_plot_jw<-unique(clust_plot_jw)
p_jw <- ggplot(clust_plot_jw, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_jw$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) +labs(title = "Jarro Winkler distance clusters")


clust_plot_cosine <- benchmark_aff_clust_cosine %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) 
clust_plot_cosine<-unique(clust_plot_cosine)
p_cosine <- ggplot(clust_plot_cosine, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_cosine$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05)+labs(title = "Cosine distance clusters")



# Global Plot
global_clust_plot_lv <- nested_affinity_cluster_lv %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
global_clust_plot_jw <- nested_affinity_cluster_jw %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
global_clust_plot_cosine <- nested_affinity_cluster_cosine %>% dplyr::select(Cluster_ID,mean_silo_score,cluster_member_count) %>% dplyr::distinct()

plt_global_lv <- ggplot(global_clust_plot_lv, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_lv$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Levenshtein distance clusters")

plt_global_jw <- ggplot(global_clust_plot_jw, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_jw$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Jarro Winkler distance clusters")

plt_global_cosine <- ggplot(global_clust_plot_cosine, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_cosine$Cluster_ID,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Cosine distance clusters")

# Save plot
ggsave(plt_global_lv, filename = paste(plots_dir,"/plt_global_lv.pdf",sep=""), height = 30, width = 21, units = "cm")
ggsave(plt_global_jw, filename = paste(plots_dir,"/plt_global_jw.pdf",sep=""), height = 30, width = 21, units = "cm")
ggsave(plt_global_cosine, filename = paste(plots_dir,"/plt_global_cosine.pdf",sep=""), height = 30, width = 21, units = "cm")


# Sillohoutte plot
benchmark_aff_clust_lv <- benchmark_aff_clust_lv[order(benchmark_aff_clust_lv$Cluster_ID,-benchmark_aff_clust_lv$silhouette_score),]
benchmark_aff_clust_jw <- benchmark_aff_clust_jw[order(benchmark_aff_clust_jw$Cluster_ID,-benchmark_aff_clust_jw$silhouette_score),]
benchmark_aff_clust_cosine <- benchmark_aff_clust_cosine[order(benchmark_aff_clust_cosine$Cluster_ID,-benchmark_aff_clust_cosine$silhouette_score),]

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_lv$Cluster_ID),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("Cluster_ID","color_density")
benchmark_aff_clust_lv<- benchmark_aff_clust_lv %>% dplyr::left_join(color_density_df,by="Cluster_ID")
benchmark_aff_clust_lv$color_density<- as.double(benchmark_aff_clust_lv$color_density)

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_jw$Cluster_ID),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("Cluster_ID","color_density")
benchmark_aff_clust_jw<- benchmark_aff_clust_jw %>% dplyr::left_join(color_density_df,by="Cluster_ID")
benchmark_aff_clust_jw$color_density<- as.double(benchmark_aff_clust_jw$color_density)

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_cosine$Cluster_ID),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("Cluster_ID","color_density")
benchmark_aff_clust_cosine<- benchmark_aff_clust_cosine %>% dplyr::left_join(color_density_df,by="Cluster_ID")
benchmark_aff_clust_cosine$color_density<- as.double(benchmark_aff_clust_cosine$color_density)

barplot(height=benchmark_aff_clust_lv$silhouette_score,names=benchmark_aff_clust_lv$Cluster_ID,ylim =c(-1,1),density=benchmark_aff_clust_lv$color_density,main="Levenshtein")
barplot(height=benchmark_aff_clust_jw$silhouette_score,names=benchmark_aff_clust_jw$Cluster_ID,ylim =c(-1,1),density=benchmark_aff_clust_jw$color_density, main="Jarro Winkler")
barplot(height=benchmark_aff_clust_cosine$silhouette_score,names=benchmark_aff_clust_cosine$Cluster_ID,ylim =c(-1,1),density=benchmark_aff_clust_cosine$color_density,main="Cosine")





save.image(file = paste(results_dir,"/script3_5thed.RData",sep=""))

write.csv(benchmark_aff_clust_lv,paste(results_dir,"/benchmark_aff_clust_lv.csv",sep=""))
write.csv(benchmark_aff_clust_jw,paste(results_dir,"/benchmark_aff_clust_jw.csv",sep=""))
write.csv(benchmark_aff_clust_cosine,paste(results_dir,"/benchmark_aff_clust_cosine.csv",sep=""))


write.csv(nested_affinity_cluster_lv,paste(results_dir,"/nested_affinity_cluster_lv.csv",sep=""))
write.csv(nested_affinity_cluster_jw,paste(results_dir,"/nested_affinity_cluster_jw.csv",sep=""))
write.csv(nested_affinity_cluster_cosine,paste(results_dir,"/nested_affinity_cluster_cosine.csv",sep=""))




