# Script computes edit distances between clinical trials tumors, WHO tumors, and NCIT tumors. 
# Performs Affinity Cluster with 3 levels of nesting. 

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
results_dir <- file.path(analysis_dir,"results")
plots_dir <- file.path(root_dir,"plots")
source(paste(util_dir,"/string_normalizing.R",sep=""))
source(paste(util_dir,"/nested_clust_edit_dist.R",sep=""))
source(paste(util_dir,"/compute_silhouette.R",sep=""))
source(paste(util_dir,"/edit_distance_nested_cluster.R",sep=""))



# Load data

load(paste(intermediate_dir,"/dissimilarity_matrix_lv.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_jw.RData",sep=""))
load(paste(intermediate_dir,"/dissimilarity_matrix_cosine.RData",sep=""))


cluster_results_lv<- read.csv(paste(results_dir,"/cluster_lv.csv",sep=""))
cluster_results_jw<- read.csv(paste(results_dir,"/cluster_jw.csv",sep=""))
cluster_results_cosine<- read.csv(paste(results_dir,"/cluster_cosine.csv",sep=""))


# Compute Simmilarity matrix for each edit distance
simmilarity_matrix_cosine = 1 - dissimilarity_matrix_cosine
simmilarity_matrix_jw = 1-dissimilarity_matrix_jw



cl <- makeCluster(25, outfile="")
registerDoParallel(cl)
df_tumor_names<- colnames(dissimilarity_matrix_lv)
normalizing_matrix_lv <- as.data.frame(matrix(nrow=length(df_tumor_names),ncol=length(df_tumor_names)))
rownames(normalizing_matrix_lv)<-df_tumor_names
colnames(normalizing_matrix_lv)<-df_tumor_names


normalizing_matrix_lv<-foreach(iter=1:length(df_tumor_names),.combine=rbind) %dopar% {
  print(iter)
  disease_name <- colnames(dissimilarity_matrix_lv)[iter]
  norm_factors<-unlist(lapply(df_tumor_names,string_normalzing,S2=disease_name))
}
rownames(normalizing_matrix_lv) <- df_tumor_names
colnames(normalizing_matrix_lv) <- df_tumor_names

simmilarity_matrix_lv <- 1- (dissimilarity_matrix_lv/normalizing_matrix_lv)
save(simmilarity_matrix_lv, file=paste(intermediate_dir,"/simmilarity_matrix_lv.RData",sep=""))

stopCluster(cl)

######### Cluster with LV ########
apclust_lv <- apcluster(simmilarity_matrix_lv) #10:26pm -12:45 pm #2:07 pm-8:00 pm
affinity_cluster_lv_df<-as.data.frame(matrix(nrow=1,ncol=2))
colnames(affinity_cluster_lv_df)<-c("Tumor_Names","Cluster_ID")
for (iter in 1: length(apclust_lv@clusters)){
  affinity_cluster_lv_df[iter,1] <- paste(names(unlist(apclust_lv@clusters[iter])),collapse = "@")
  affinity_cluster_lv_df[iter,2] <- iter
}
affinity_cluster_lv_df<- affinity_cluster_lv_df %>% separate_rows(Tumor_Names, sep = '@')

######### Cluster with jw ########
apclust_jw <- apcluster(simmilarity_matrix_jw) #8:32 pm -10:06 pm. 8:00 pm 
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
nested_affinity_cluster_lv2<- edit_distance_nested_cluster(affinity_cluster_lv_df,simmilarity_matrix_lv)
nested_affinity_cluster_jw2<- edit_distance_nested_cluster(affinity_cluster_jw_df,simmilarity_matrix_jw)
nested_affinity_cluster_cosine2<- edit_distance_nested_cluster(affinity_cluster_cosine_df,simmilarity_matrix_cosine)


### Nested LV Clustering

nested_affinity_cluster_lv <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_lv_df,
                                                     dist_mat = simmilarity_matrix_lv)

### Nested JW Clustering

nested_affinity_cluster_jw <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_jw_df,
                                                     dist_mat = simmilarity_matrix_jw)
### Nested Cosine Clustering

nested_affinity_cluster_cosine <- nested_clust_edit_dist(n=3,affinity_cluster_df = affinity_cluster_cosine_df,
                                                     dist_mat = simmilarity_matrix_cosine)
# Compute Silo Dist
#nested_affinity_cluster_lv<-compute_silhouette(cluster_df = nested_affinity_cluster_lv,dist_mat = dissimilarity_matrix_lv)
#nested_affinity_cluster_jw<-compute_silhouette(cluster_df = nested_affinity_cluster_jw,dist_mat = dissimilarity_matrix_jw)
#nested_affinity_cluster_cosine<-compute_silhouette(cluster_df = nested_affinity_cluster_cosine,dist_mat = dissimilarity_matrix_cosine)

nested_affinity_cluster_lv2<-compute_silhouette(cluster_df = nested_affinity_cluster_lv2,dist_mat = dissimilarity_matrix_lv)
nested_affinity_cluster_jw2<-compute_silhouette(cluster_df = nested_affinity_cluster_jw2,dist_mat = dissimilarity_matrix_jw)
nested_affinity_cluster_cosine2<-compute_silhouette(cluster_df = nested_affinity_cluster_cosine2,dist_mat = dissimilarity_matrix_cosine)


mean_freq_lv <- nested_affinity_cluster_lv %>% dplyr::select(SubsetCluster_IDs, silhouette_score) %>% dplyr::group_by(SubsetCluster_IDs) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n()) 
mean_freq_jw <- nested_affinity_cluster_jw %>% dplyr::select(SubsetCluster_IDs, silhouette_score) %>% dplyr::group_by(SubsetCluster_IDs) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n())
mean_freq_cosine <- nested_affinity_cluster_cosine %>% dplyr::select(SubsetCluster_IDs, silhouette_score) %>% dplyr::group_by(SubsetCluster_IDs) %>% dplyr::summarise(mean_silo_score=mean(silhouette_score),cluster_member_count =dplyr::n())


nested_affinity_cluster_lv<- nested_affinity_cluster_lv %>% dplyr::left_join(mean_freq_lv,by="SubsetCluster_IDs")
nested_affinity_cluster_jw<- nested_affinity_cluster_jw %>% dplyr::left_join(mean_freq_jw,by="SubsetCluster_IDs")
nested_affinity_cluster_cosine<- nested_affinity_cluster_cosine %>% dplyr::left_join(mean_freq_cosine,by="SubsetCluster_IDs")

# Select benchmarks 
benchmark_tumors <- c("b cell lymphoma", "neuroblastoma", "triple negative breast cancer",
                      "unresectable lung carcinoma", "liposarcoma","cancer of the liver",
                      "smoldering myeloma")

subcluster_lv <- nested_affinity_cluster_lv$SubsetCluster_IDs[nested_affinity_cluster_lv$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_lv <- nested_affinity_cluster_lv %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_lv)

subcluster_jw <- nested_affinity_cluster_jw$SubsetCluster_IDs[nested_affinity_cluster_jw$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_jw <- nested_affinity_cluster_jw %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_jw)

subcluster_cosine <- nested_affinity_cluster_cosine$SubsetCluster_IDs[nested_affinity_cluster_cosine$Tumor_Names %in% benchmark_tumors]
benchmark_aff_clust_cosine <- nested_affinity_cluster_cosine %>% dplyr::filter(SubsetCluster_IDs %in% subcluster_cosine)


# Plot 

clust_plot_lv <- benchmark_aff_clust_lv %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) 
clust_plot_lv<-unique(clust_plot_lv)
p_lv <- ggplot(clust_plot_lv, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_lv$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Levenshtein distance clusters")


clust_plot_jw <- benchmark_aff_clust_jw %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) 
clust_plot_jw<-unique(clust_plot_jw)
p_jw <- ggplot(clust_plot_jw, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_jw$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) +labs(title = "Jarro Winkler distance clusters")


clust_plot_cosine <- benchmark_aff_clust_cosine %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) 
clust_plot_cosine<-unique(clust_plot_cosine)
p_cosine <- ggplot(clust_plot_cosine, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=clust_plot_cosine$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05)+labs(title = "Cosine distance clusters")



# Global Plot
global_clust_plot_lv <- nested_affinity_cluster_lv %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
global_clust_plot_jw <- nested_affinity_cluster_jw %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) %>% dplyr::distinct()
global_clust_plot_cosine <- nested_affinity_cluster_cosine %>% dplyr::select(SubsetCluster_IDs,mean_silo_score,cluster_member_count) %>% dplyr::distinct()

plt_global_lv <- ggplot(global_clust_plot_lv, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_lv$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Levenshtein distance clusters")

plt_global_jw <- ggplot(global_clust_plot_jw, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_jw$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Jarro Winkler distance clusters")

plt_global_cosine <- ggplot(global_clust_plot_cosine, aes(x=cluster_member_count, y=mean_silo_score)) +geom_point() +
  geom_text(label=global_clust_plot_cosine$SubsetCluster_IDs,check_overlap = TRUE,angle = 45,vjust = 0, nudge_y = 0.05) + labs(title = "Cosine distance clusters")

# Save plot
ggsave(plt_global_lv, filename = paste(plots_dir,"/plt_global_lv.pdf",sep=""), height = 30, width = 21, units = "cm")
ggsave(plt_global_jw, filename = paste(plots_dir,"/plt_global_jw.pdf",sep=""), height = 30, width = 21, units = "cm")
ggsave(plt_global_cosine, filename = paste(plots_dir,"/plt_global_cosine.pdf",sep=""), height = 30, width = 21, units = "cm")


# Sillohoutte plot
benchmark_aff_clust_lv <- benchmark_aff_clust_lv[order(benchmark_aff_clust_lv$SubsetCluster_IDs,-benchmark_aff_clust_lv$silhouette_score),]
benchmark_aff_clust_jw <- benchmark_aff_clust_jw[order(benchmark_aff_clust_jw$SubsetCluster_IDs,-benchmark_aff_clust_jw$silhouette_score),]
benchmark_aff_clust_cosine <- benchmark_aff_clust_cosine[order(benchmark_aff_clust_cosine$SubsetCluster_IDs,-benchmark_aff_clust_cosine$silhouette_score),]

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_lv$SubsetCluster_IDs),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("SubsetCluster_IDs","color_density")
benchmark_aff_clust_lv<- benchmark_aff_clust_lv %>% dplyr::left_join(color_density_df,by="SubsetCluster_IDs")
benchmark_aff_clust_lv$color_density<- as.double(benchmark_aff_clust_lv$color_density)

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_jw$SubsetCluster_IDs),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("SubsetCluster_IDs","color_density")
benchmark_aff_clust_jw<- benchmark_aff_clust_jw %>% dplyr::left_join(color_density_df,by="SubsetCluster_IDs")
benchmark_aff_clust_jw$color_density<- as.double(benchmark_aff_clust_jw$color_density)

color_density_df <- as.data.frame(cbind(unique(benchmark_aff_clust_cosine$SubsetCluster_IDs),c(1,10,30,50,70,90,120)))
colnames(color_density_df)<- c("SubsetCluster_IDs","color_density")
benchmark_aff_clust_cosine<- benchmark_aff_clust_cosine %>% dplyr::left_join(color_density_df,by="SubsetCluster_IDs")
benchmark_aff_clust_cosine$color_density<- as.double(benchmark_aff_clust_cosine$color_density)

barplot(height=benchmark_aff_clust_lv$silhouette_score,names=benchmark_aff_clust_lv$SubsetCluster_IDs,ylim =c(-1,1),density=benchmark_aff_clust_lv$color_density,main="Levenshtein")
barplot(height=benchmark_aff_clust_jw$silhouette_score,names=benchmark_aff_clust_jw$SubsetCluster_IDs,ylim =c(-1,1),density=benchmark_aff_clust_jw$color_density, main="Jarro Winkler")
barplot(height=benchmark_aff_clust_cosine$silhouette_score,names=benchmark_aff_clust_cosine$SubsetCluster_IDs,ylim =c(-1,1),density=benchmark_aff_clust_cosine$color_density,main="Cosine")


## Write workspace
save.image(file = "script3.RData")

write.csv(benchmark_aff_clust_lv,paste(results_dir,"/benchmark_aff_clust_lv.csv",sep=""))
write.csv(benchmark_aff_clust_jw,paste(results_dir,"/benchmark_aff_clust_jw.csv",sep=""))
write.csv(benchmark_aff_clust_cosine,paste(results_dir,"/benchmark_aff_clust_cosine.csv",sep=""))


write.csv(nested_affinity_cluster_lv,paste(results_dir,"/nested_affinity_cluster_lv.csv",sep=""))
write.csv(nested_affinity_cluster_jw,paste(results_dir,"/nested_affinity_cluster_jw.csv",sep=""))
write.csv(nested_affinity_cluster_cosine,paste(results_dir,"/nested_affinity_cluster_cosine.csv",sep=""))

write.csv(dissimilarity_matrix_cosine,paste(intermediate_dir,"/dissimilarity_matrix_cosine.csv",sep=""))
write.csv(dissimilarity_matrix_jw,paste(intermediate_dir,"/dissimilarity_matrix_jw.csv",sep=""))
write.csv(dissimilarity_matrix_lv,paste(intermediate_dir,"/dissimilarity_matrix_lv.csv",sep=""))