#This script is used to generate the box plots for figure 1 and 2 . The box-plots compare the standardization results for the methods LTE-3 + Euclidean Dist and all-MiniLM-L12-v2+Euclidean Dist.  

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

source(paste(util_dir,"/find_euclidean_match.R",sep=""))
# Load the annotations for 5th edition 
#tumor_5th_edition<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))
#tumor_all_edition<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))

#tumor_all_edition<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
#tumor_5th_edition<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_edition<-read_excel(paste(result_dir,"/tumor_manually_validated_all_corrected_May23.xlsx",sep = ""))
tumor_5th_edition<-read_excel(paste(result_dir_5th,"/tumor_manually_validated_5th_corrected_May23.xlsx",sep = ""))





# Pick only the Euclidean distance standardization with V3 embeddings
tumor_5th_edition<-tumor_5th_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3,
                                                       euclidean_dist_MiniLM_L12_v2,valid_euclidean_dist_MiniLM_L12_v2,
)
tumor_all_edition<-tumor_all_edition %>% dplyr::select(nct_id,Tumor_Names,ground_truth_val,ground_truth,
                                                       euclidean_dist_v3,valid_euclidean_dist_v3,
                                                       euclidean_dist_MiniLM_L12_v2,valid_euclidean_dist_MiniLM_L12_v2,
)


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
MiniLM_L12_v2_embeddings<-MiniLM_L12_v2_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")


WHO_5th_edition<-find_euclidean_match(WHO_5th_edition,MiniLM_L12_v2_embeddings,tumor_5th_edition[,c("Tumor_Names","euclidean_dist_MiniLM_L12_v2")],c("MiniLM_L12_v2_Matches","Euclidean_Dist_MiniLM_L12_v2"))
WHO_all_edition<-find_euclidean_match(WHO_all_edition,MiniLM_L12_v2_embeddings,tumor_all_edition[,c("Tumor_Names","euclidean_dist_MiniLM_L12_v2")],c("MiniLM_L12_v2_Matches","Euclidean_Dist_MiniLM_L12_v2"))









# Join 5th edition data
tumor_5th_edition<- tumor_5th_edition %>% dplyr::left_join(WHO_5th_edition%>%dplyr::select(Tumor_Names,Euclidean_Dist_LTE3,Euclidean_Dist_MiniLM_L12_v2),by="Tumor_Names")

#Join all edition data
tumor_all_edition<- tumor_all_edition %>% dplyr::left_join(WHO_all_edition%>%dplyr::select(Tumor_Names,Euclidean_Dist_LTE3,Euclidean_Dist_MiniLM_L12_v2),by="Tumor_Names")




# Remove cases where there were no ground truths (NF = Not Found)
tumor_5th_edition<- tumor_5th_edition %>%filter(ground_truth !="NF")
tumor_all_edition<- tumor_all_edition %>%filter(ground_truth !="NF")


# data table for 5th and all edition Euclidean V3 distance
distances_5th_edition<- tumor_5th_edition%>%dplyr::select(Euclidean_Dist_LTE3,valid_euclidean_dist_v3, Euclidean_Dist_MiniLM_L12_v2, valid_euclidean_dist_MiniLM_L12_v2)
distances_all_edition<- tumor_all_edition%>%dplyr::select(Euclidean_Dist_LTE3,valid_euclidean_dist_v3, Euclidean_Dist_MiniLM_L12_v2, valid_euclidean_dist_MiniLM_L12_v2)

distances_5th_edition<-distances_5th_edition %>% mutate(standardization_result_LTE3= case_when(valid_euclidean_dist_v3==1~ "Correctly Standardized",
                                                                                               valid_euclidean_dist_v3==0~"Incorrectly Standardized"))

distances_5th_edition<-distances_5th_edition %>% mutate(standardization_result_all_MiniLM_L12_v2= case_when(valid_euclidean_dist_MiniLM_L12_v2==1~ "Correctly Standardized",
                                                                                                            valid_euclidean_dist_MiniLM_L12_v2==0~"Incorrectly Standardized"))




distances_all_edition<-distances_all_edition %>% mutate(standardization_result_LTE3= case_when(valid_euclidean_dist_v3==1~ "Correctly Standardized",
                                                                                               valid_euclidean_dist_v3==0~"Incorrectly Standardized"))

distances_all_edition<-distances_all_edition %>% mutate(standardization_result_all_MiniLM_L12_v2= case_when(valid_euclidean_dist_MiniLM_L12_v2==1~ "Correctly Standardized",
                                                                                                            valid_euclidean_dist_MiniLM_L12_v2==0~"Incorrectly Standardized"))


colnames(distances_5th_edition)[c(1,3)]<-c("Euclidean_distance_LTE_embedding","Euclidean_distance_all_MiniLM_L12_v2_embedding")
colnames(distances_all_edition)[c(1,3)]<-c("Euclidean_distance_LTE_embedding","Euclidean_distance_all_MiniLM_L12_v2_embedding")


Plt_5th_ed_LTE3<- ggplot(distances_5th_edition, aes(x=standardization_result_LTE3, y=Euclidean_distance_LTE_embedding,
                                                    color=standardization_result_LTE3)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ 
  labs(y= "Euclidean distance in LTE-3 embedding space", x = "Standardization results for WHO 5th edition,n=1044",color = "Standardization Result LTE-3+Euclidean distance")+ggtitle("a")+
  theme(plot.title = element_text(size = 16, face = "bold"),axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),legend.title = element_text(size = 16),axis.title.x = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16))

Plt_all_ed_LTE3<- ggplot(distances_all_edition, aes(x=standardization_result_LTE3, y=Euclidean_distance_LTE_embedding,
                                                    color=standardization_result_LTE3)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ 
  labs(y= "Euclidean distance in LTE-3 embedding space", x = "Standardization results for WHO all editions,n=1118",color = "Standardization Result LTE-3+Euclidean distance")+ggtitle("b")+
  theme(plot.title = element_text(size = 16, face = "bold"),axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),legend.title = element_text(size = 16),axis.title.x = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16))



Plt_5th_ed_MiniLM_L12_v2<- ggplot(distances_5th_edition, aes(x=standardization_result_all_MiniLM_L12_v2, y=Euclidean_distance_all_MiniLM_L12_v2_embedding,
                                                             color=standardization_result_all_MiniLM_L12_v2)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+
  labs(y= "Euclidean distance in all-MiniLM-L12-v2 embedding space", x = "Standardization results for WHO 5th edition,n=1044",color = "Standardization Result all-MiniLM-L12-v2+Euclidean distance")+ggtitle("a")+
  theme(plot.title = element_text(size = 16, face = "bold"),axis.title = element_text(size = 16),axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

Plt_all_ed_MiniLM_L12_v2<- ggplot(distances_all_edition, aes(x=standardization_result_all_MiniLM_L12_v2, y=Euclidean_distance_all_MiniLM_L12_v2_embedding,
                                                             color=standardization_result_all_MiniLM_L12_v2)) + geom_boxplot()+ scale_fill_brewer(palette="Dark2")+ 
  labs(y= "Euclidean distance in all-MiniLM-L12-v2 embedding space", x = "Standardization results for WHO all editions,n=1118",color = "Standardization Result all-MiniLM-L12-v2+Euclidean distance")+ggtitle("b")+
  theme(plot.title = element_text(size = 16, face = "bold"),axis.title = element_text(size = 16),axis.title.x = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 16),legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))




distances_5th_edition_LTE_correct<- tumor_5th_edition%>%filter(valid_euclidean_dist_v3==1)
distances_5th_edition_LTE_wrong<- tumor_5th_edition%>%filter(valid_euclidean_dist_v3==0)
distances_all_edition_LTE_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==1)
distances_all_edition_LTE_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_v3==0)

distances_5th_edition_MiniLM_L12_v2_correct<- tumor_5th_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==1)
distances_5th_edition_MiniLM_L12_v2_wrong<- tumor_5th_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==0)
distances_all_edition_MiniLM_L12_v2_correct<- tumor_all_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==1)
distances_all_edition_MiniLM_L12_v2_wrong<- tumor_all_edition%>%filter(valid_euclidean_dist_MiniLM_L12_v2==0)



summary_5th_LTE_correct<-signif(summary(distances_5th_edition_LTE_correct$Euclidean_Dist_LTE3),digits=3)
summary_5th_LTE_wrong<-signif(summary(distances_5th_edition_LTE_wrong$Euclidean_Dist_LTE3),digits=3)

summary_all_LTE_correct<-signif(summary(distances_all_edition_LTE_correct$Euclidean_Dist_LTE3),digits=3)
summary_all_LTE_wrong<-signif(summary(distances_all_edition_LTE_wrong$Euclidean_Dist_LTE3),digits=3)



summary_5th_MiniLM_L12_v2_correct<-signif(summary(distances_5th_edition_MiniLM_L12_v2_correct$Euclidean_Dist_MiniLM_L12_v2),digits=3)
summary_5th_MiniLM_L12_v2_wrong<-signif(summary(distances_5th_edition_MiniLM_L12_v2_wrong$Euclidean_Dist_MiniLM_L12_v2),digits=3)

summary_all_MiniLM_L12_v2_correct<-signif(summary(distances_all_edition_MiniLM_L12_v2_correct$Euclidean_Dist_MiniLM_L12_v2),digits=3)
summary_all_MiniLM_L12_v2_wrong<-signif(summary(distances_all_edition_MiniLM_L12_v2_wrong$Euclidean_Dist_MiniLM_L12_v2),digits=3)







print("LTE-3 Statistic for all editions")
print(summary_all_LTE_correct)
print(summary_all_LTE_wrong)

print("LTE-3 Statistic for 5th editions")
print(summary_5th_LTE_correct)
print(summary_5th_LTE_wrong)

print("MiniLM_L12_v2 Statistic for all editions")
print(summary_all_MiniLM_L12_v2_correct)
print(summary_all_MiniLM_L12_v2_wrong)

print("MiniLM_L12_v2 Statistic for 5th editions")
print(summary_5th_MiniLM_L12_v2_correct)
print(summary_5th_MiniLM_L12_v2_wrong)




Plt_LTE<-ggarrange(Plt_5th_ed_LTE3, Plt_all_ed_LTE3,nrow = 1,ncol = 2)
ggsave(paste(root_dir,"/Paper/MLWA/High_Res_Fig/figure5_revise.jpg",sep=""), dpi = 300, width = 70, height = 25, units = "cm")

Plt_MiniLM_L12_v2<-ggarrange(Plt_5th_ed_MiniLM_L12_v2, Plt_all_ed_MiniLM_L12_v2,nrow = 1,ncol = 2)
ggsave(paste(root_dir,"/Paper/MLWA/High_Res_Fig/figure6_revise.jpg",sep=""), dpi = 300, width = 70, height = 25, units = "cm")
