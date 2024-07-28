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

save.image(file = "script3_5thed.RData")




