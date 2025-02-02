#This script is used to generate the box plots for figure 1 . The box-plots compare the standardization results of LTE-3 + Euclidean Dist method.

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(ghql)
  library(readxl)
  library(dbscan)
  library(isotree)
  library(ggpubr)
  
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
intermediate_dir_5th <- file.path(analysis_dir,"intermediate_5th")
result_dir <- file.path(analysis_dir,"results")
result_dir_5th <- file.path(analysis_dir,"results_5th")


# Load the annotations for 5th edition 
#tumor_5th_edition<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))
#tumor_all_edition<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))

tumor_all_edition<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5th_edition<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

# Pick only the Euclidean distance standardization with V3 embeddings
tumor_5th_edition<-tumor_5th_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3,
                                                       euclidean_dist_MiniLM_L12_v2,valid_euclidean_dist_MiniLM_L12_v2,
                                                       euclidean_dist_e5_large,valid_euclidean_dist_e5_large)
tumor_all_edition<-tumor_all_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3,
                                                       euclidean_dist_MiniLM_L12_v2,valid_euclidean_dist_MiniLM_L12_v2,
                                                       euclidean_dist_e5_large,valid_euclidean_dist_e5_large)


# Load 5th edition distance for LTE-3
affinity_cluster_v3_reassigned_5thed_df<-read.csv(paste(intermediate_dir_5th,"/affinity_cluster_v3_reassigned_df_5thed.csv",sep=""))
WHO_5th_edition<-affinity_cluster_v3_reassigned_5thed_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance)
colnames(WHO_5th_edition)[c(2,3)]<- c("LTE3_Matches","Euclidean_Dist_LTE3")

# Load all edition distance for LTE-3
affinity_cluster_v3_reassigned_df<-read.csv(paste(intermediate_dir,"/affinity_cluster_v3_reassigned_df.csv",sep=""))
WHO_all_edition<-affinity_cluster_v3_reassigned_df %>% dplyr::select(Tumor_Names,WHO_Matches,WHO_distance)
colnames(WHO_all_edition)[c(2,3)]<- c("LTE3_Matches","Euclidean_Dist_LTE3")

WHO_5th_edition<-WHO_5th_edition%>%filter(Tumor_Names %in% tumor_5th_edition$Tumor_Names)
WHO_all_edition<-WHO_all_edition%>%filter(Tumor_Names %in% tumor_all_edition$Tumor_Names)


# Load the embeddings for MiniLM_V12 and E-5_Large
MiniLM_L12_v2_embeddings<-read.csv(paste(data_dir,"/Embeddings/all-MiniLM-L12-v2.csv",sep=""))
e5_large_embeddings<-read.csv(paste(data_dir,"/Embeddings/e5-large.csv",sep="")) 
MiniLM_L12_v2_embeddings<-MiniLM_L12_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
e5_large_embeddings<-e5_large_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")





WHO_5th_edition<-find_euclidean_match(WHO_5th_edition,MiniLM_L12_v2_embeddings,tumor_5th_edition[,c("Tumor_Names","euclidean_dist_MiniLM_L12_v2")],c("MiniLM_L12_v2_Matches","Euclidean_Dist_MiniLM_L12_v2"))
WHO_5th_edition<-find_euclidean_match(WHO_5th_edition,e5_large_embeddings,tumor_5th_edition[,c("Tumor_Names","euclidean_dist_e5_large")],c("e5large_Matches","Euclidean_Dist_e5_large"))




WHO_all_edition<-find_euclidean_match(WHO_all_edition,MiniLM_L12_v2_embeddings,tumor_all_edition[,c("Tumor_Names","euclidean_dist_MiniLM_L12_v2")],c("MiniLM_L12_v2_Matches","Euclidean_Dist_MiniLM_L12_v2"))
WHO_all_edition<-find_euclidean_match(WHO_all_edition,e5_large_embeddings,tumor_all_edition[,c("Tumor_Names","euclidean_dist_e5_large")],c("e5large_Matches","Euclidean_Dist_e5_large"))









# Join 5th edition data
tumor_5th_edition<- tumor_5th_edition %>% dplyr::left_join(WHO_5th_edition%>%dplyr::select(Tumor_Names,Euclidean_Dist_LTE3,Euclidean_Dist_MiniLM_L12_v2,Euclidean_Dist_e5_large),by="Tumor_Names")

#Join all edition data
tumor_all_edition<- tumor_all_edition %>% dplyr::left_join(WHO_all_edition%>%dplyr::select(Tumor_Names,Euclidean_Dist_LTE3,Euclidean_Dist_MiniLM_L12_v2,Euclidean_Dist_e5_large),by="Tumor_Names")




# Remove cases where there were no ground truths (NF = Not Found)
tumor_5th_edition<- tumor_5th_edition %>%filter(ground_truth !="NF")
tumor_all_edition<- tumor_all_edition %>%filter(ground_truth !="NF")


# data table for 5th and all edition Euclidean V3 distance
distances_5th_edition<- tumor_5th_edition%>%dplyr::select(Euclidean_Dist_LTE3,valid_euclidean_dist_v3, Euclidean_Dist_MiniLM_L12_v2, valid_euclidean_dist_MiniLM_L12_v2, Euclidean_Dist_e5_large, valid_euclidean_dist_e5_large)
distances_all_edition<- tumor_all_edition%>%dplyr::select(Euclidean_Dist_LTE3,valid_euclidean_dist_v3, Euclidean_Dist_MiniLM_L12_v2, valid_euclidean_dist_MiniLM_L12_v2, Euclidean_Dist_e5_large, valid_euclidean_dist_e5_large)

distances_5th_edition<-distances_5th_edition %>% mutate(standardization_result_LTE3= case_when(valid_euclidean_dist_v3==1~ "Correctly Standardized",
                                                                                         valid_euclidean_dist_v3==0~"Incorrectly Standardized"))

distances_5th_edition<-distances_5th_edition %>% mutate(standardization_result_MiniLM_L12_v2= case_when(valid_euclidean_dist_MiniLM_L12_v2==1~ "Correctly Standardized",
                                                                                               valid_euclidean_dist_MiniLM_L12_v2==0~"Incorrectly Standardized"))

distances_5th_edition<-distances_5th_edition %>% mutate(standardization_result_e5_large= case_when(valid_euclidean_dist_e5_large==1~ "Correctly Standardized",
                                                                                                        valid_euclidean_dist_e5_large==0~"Incorrectly Standardized"))

                                                        


distances_all_edition<-distances_all_edition %>% mutate(standardization_result_LTE3= case_when(valid_euclidean_dist_v3==1~ "Correctly Standardized",
                                                                                               valid_euclidean_dist_v3==0~"Incorrectly Standardized"))

distances_all_edition<-distances_all_edition %>% mutate(standardization_result_MiniLM_L12_v2= case_when(valid_euclidean_dist_MiniLM_L12_v2==1~ "Correctly Standardized",
                                                                                                        valid_euclidean_dist_MiniLM_L12_v2==0~"Incorrectly Standardized"))

distances_all_edition<-distances_all_edition %>% mutate(standardization_result_e5_large= case_when(valid_euclidean_dist_e5_large==1~ "Correctly Standardized",
                                                                                                   valid_euclidean_dist_e5_large==0~"Incorrectly Standardized"))

colnames(distances_5th_edition)[c(1,3,5)]<-c("Euclidean_distance_LTE_embedding","Euclidean_distance_MiniLM_L12_v2_embedding","Euclidean_distance_E5_Large_embedding")
colnames(distances_all_edition)[c(1,3,5)]<-c("Euclidean_distance_LTE_embedding","Euclidean_distance_MiniLM_L12_v2_embedding","Euclidean_distance_E5_Large_embedding")


Plt_5th_ed_LTE3<- ggplot(distances_5th_edition, aes(x=standardization_result_LTE3, y=Euclidean_distance_LTE_embedding,
                                  color=standardization_result_LTE3)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ labs(y= "Euclidean distance in LTE-3 embedding space", x = "Standardization results for WHO 5th edition")+ggtitle("a")

Plt_all_ed_LTE3<- ggplot(distances_all_edition, aes(x=standardization_result_LTE3, y=Euclidean_distance_LTE_embedding,
                                              color=standardization_result_LTE3)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2") + labs(y= "Euclidean distance in LTE-3 embedding space", x = "Standardization results for WHO all edition")+ggtitle("b")



Plt_5th_ed_MiniLM_L12_v2<- ggplot(distances_5th_edition, aes(x=standardization_result_MiniLM_L12_v2, y=Euclidean_distance_MiniLM_L12_v2_embedding,
                                                    color=standardization_result_MiniLM_L12_v2)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ labs(y= "Euclidean distance in MiniLM_L12_v2 embedding space", x = "Standardization results for WHO 5th edition")+ggtitle("a")

Plt_all_ed_MiniLM_L12_v2<- ggplot(distances_all_edition, aes(x=standardization_result_MiniLM_L12_v2, y=Euclidean_distance_MiniLM_L12_v2_embedding,
                                                    color=standardization_result_MiniLM_L12_v2)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2") + labs(y= "Euclidean distance in MiniLM_L12_v2 embedding space", x = "Standardization results for WHO all edition")+ggtitle("b")




Plt_5th_ed_e5_large<- ggplot(distances_5th_edition, aes(x=standardization_result_e5_large, y=Euclidean_distance_E5_Large_embedding,
                                                    color=standardization_result_e5_large)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ labs(y= "Euclidean distance in E5_Large embedding space", x = "Standardization results for WHO 5th edition")+ggtitle("a")

Plt_all_ed_e5_large<- ggplot(distances_all_edition, aes(x=standardization_result_e5_large, y=Euclidean_distance_E5_Large_embedding,
                                                    color=standardization_result_e5_large)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2") + labs(y= "Euclidean distance in E5_Large embedding space", x = "Standardization results for WHO all edition")+ggtitle("b")



distances_5th_edition_LTE_correct<- tumor_5th_edition%>%filter(valid_euclidean_dist_v3==1)
distances_5th_edition_LTE_wrong<- tumor_5th_edition%>%filter(valid_euclidean_dist_v3==0)
distances_all_edition_LTE_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==1)
distances_all_edition_LTE_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==0)

distances_5th_edition_MiniLM_L12_v2_correct<- tumor_5th_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==1)
distances_5th_edition_MiniLM_L12_v2_wrong<- tumor_5th_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==0)
distances_all_edition_MiniLM_L12_v2_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==1)
distances_all_edition_MiniLM_L12_v2_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==0)

distances_5th_edition_e5_large_correct<- tumor_5th_edition%>%filter(valid_euclidean_dist_e5_large==1)
distances_5th_edition_e5_large_wrong<- tumor_5th_edition%>%filter(valid_euclidean_dist_e5_large==0)
distances_all_edition_e5_large_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_e5_large==1)
distances_all_edition_e5_large_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_e5_large==0)


summary_5th_correct <- lapply(distances_5th_edition_correct, summary)
summary_5th_wrong <- lapply(distances_5th_edition_wrong, summary)

print(summary_5th_correct$WHO_distance)
print(summary_5th_wrong$WHO_distance)



distances_all_edition_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==1)
distances_all_edition_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==0)

summary_all_correct <- lapply(distances_all_edition_correct, summary)
summary_all_wrong <- lapply(distances_all_edition_wrong, summary)

print(summary_all_correct$WHO_distance)
print(summary_all_wrong$WHO_distance)


p3<-ggarrange(Plt_5th_ed, Plt_all_ed,nrow = 1,ncol = 2)
