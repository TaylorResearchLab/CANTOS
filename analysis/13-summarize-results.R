#This script is used evaluate the accuracy of each clustering+standardization methods. 

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


#tumor_all_gt<-read.csv(paste(result_dir,"/tumor_sample_df_gt_annotated_all.csv",sep = ""))
#tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_sample_df_gt_annotated_5th.csv",sep = ""))

tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

all<-(colSums(tumor_all_gt[,c(seq(6,ncol(tumor_all_gt),2))]))/nrow(tumor_all_gt)
fifth<-(colSums(tumor_5thed_gt[,c(seq(6,ncol(tumor_5thed_gt),2))]))/nrow(tumor_5thed_gt)
all<-as.data.frame(all)
fifth<-as.data.frame(fifth)

all$method<-rownames(all)
fifth$method<-rownames(fifth)
rownames(all)<-NULL
rownames(fifth)<-NULL
all <- all %>%dplyr::select(method,all)
fifth<- fifth %>%dplyr::select(method,fifth)

all<-all[order(all$all,decreasing = TRUE),]
fifth<-fifth[order(fifth$fifth,decreasing = TRUE),]

for(iter in 1:nrow(all)){
  if(all$method[iter]=="valid_euclidean_dist_v3" ){
    all$method[iter]="LTE-3+Euclidean Dist"
  }else if(all$method[iter]=="valid_af_v3"){
    all$method[iter]="LTE-3+AP"
    
  }else if(all$method[iter]=="valid_euclidean_dist_ada2"){
    all$method[iter]="ADA-002+Euclidean Dist"
    
  }else if(all$method[iter]=="valid_af_ad2"){
    all$method[iter]="ADA-002+AP"
    
  }else if(all$method[iter]=="valid_kmeans_v3"){
    all$method[iter]="LTE-3+KMeans"
    
  }else if(all$method[iter]=="valid_kmeans_ada2"){
    all$method[iter]="ADA-002+KMeans"
    
  }else if(all$method[iter]=="valid_lv_match"){
    all$method[iter]="Levenshtein"
    
  }else if(all$method[iter]=="valid_af_lv"){
    all$method[iter]="Levenshtein+AP"
    
  }else if(all$method[iter]=="valid_jw_match"){
    all$method[iter]="Jarro Winkler"
  }else if(all$method[iter]=="valid_af_jw"){
    all$method[iter]="Jarro Winkler+AP"
    
  }else if(all$method[iter]=="valid_cosine_match"){
    all$method[iter]="Cosine"
    
  }else if(all$method[iter]=="valid_af_cosine"){
    all$method[iter]="Cosine+AP"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_llama"){
    all$method[iter]="LLama3.0+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_biobert"){
    all$method[iter]="BioBERT+Euclidean Dist"
    
  }else if(all$method[iter]=="valid_euclidean_dist_medllama"){
    all$method[iter]="MedLlama_13B+Euclidean Dist"
    
  }else if(all$method[iter]=="valid_euclidean_dist_pubmedbert"){
    all$method[iter]="PubMedBERT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_modernbert"){
    all$method[iter]="ModernBERT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_medllama_7b"){
    all$method[iter]="MedLlama2+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_llama_32_3b"){
    all$method[iter]="LLama3.2_3B+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_phi4"){
    all$method[iter]="Phi-4+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_llama_33_70b"){
    all$method[iter]="Llama3.3_70B+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_MiniLM_L6_v2"){
    all$method[iter]="all-MiniLM-L6-v2+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_mpnet_base"){
    all$method[iter]="all-mpnet-base-v2+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_e5_large"){
    all$method[iter]="e5-large+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_gtr_t5_large"){
    all$method[iter]="gtr-t5-large+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_roberta"){
    all$method[iter]="all-roberta-large-v1+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_MiniLM_L12_v2"){
    all$method[iter]="all-MiniLM-L12-v2+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_labse"){
    all$method[iter]="LaBSE+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_sciBERT"){
    all$method[iter]="sciBERT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_sapBERT"){
    all$method[iter]="sapBERT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_cohere"){
    all$method[iter]="embed-english-v2.0+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_deepseek"){
    all$method[iter]="DeepSeek_8B+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_BioGPT"){
    all$method[iter]="BioGPT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_clincalBERT"){
    all$method[iter]="ClinicalBERT+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_e5large_v2"){
    all$method[iter]="e5large_v2+Euclidean Dist"
    
  }
  else if(all$method[iter]=="valid_euclidean_dist_nomic"){
    all$method[iter]="nomic+Euclidean Dist"
    
  }
}

for(iter in 1:nrow(fifth)){
  if(fifth$method[iter]=="valid_euclidean_dist_v3" ){
    fifth$method[iter]="LTE-3+Euclidean Dist"
  }else if(fifth$method[iter]=="valid_af_v3"){
    fifth$method[iter]="LTE-3+AP"
    
  }else if(fifth$method[iter]=="valid_euclidean_dist_ada2"){
    fifth$method[iter]="ADA-002+Euclidean Dist"
    
  }else if(fifth$method[iter]=="valid_af_ad2"){
    fifth$method[iter]="ADA-002+AP"
    
  }else if(fifth$method[iter]=="valid_kmeans_v3"){
    fifth$method[iter]="LTE-3+KMeans"
    
  }else if(fifth$method[iter]=="valid_kmeans_ada2"){
    fifth$method[iter]="ADA-002+KMeans"
    
  }else if(fifth$method[iter]=="valid_lv_match"){
    fifth$method[iter]="Levenshtein"
    
  }else if(fifth$method[iter]=="valid_af_lv"){
    fifth$method[iter]="Levenshtein+AP"
    
  }else if(fifth$method[iter]=="valid_jw_match"){
    fifth$method[iter]="Jarro Winkler"
  }else if(fifth$method[iter]=="valid_af_jw"){
    fifth$method[iter]="Jarro Winkler+AP"
    
  }else if(fifth$method[iter]=="valid_cosine_match"){
    fifth$method[iter]="Cosine"
    
  }else if(fifth$method[iter]=="valid_af_cosine"){
    fifth$method[iter]="Cosine+AP"
    
  }else if(fifth$method[iter]=="valid_euclidean_dist_llama"){
    fifth$method[iter]="LLama3.0+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_biobert"){
    fifth$method[iter]="BioBERT+Euclidean Dist"
    
  }else if(fifth$method[iter]=="valid_euclidean_dist_medllama"){
    fifth$method[iter]="MedLlama_13b+Euclidean Dist"
    
  }else if(fifth$method[iter]=="valid_euclidean_dist_pubmedbert"){
    fifth$method[iter]="PubMedBERT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_modernbert"){
    fifth$method[iter]="ModernBERT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_medllama_7b"){
    fifth$method[iter]="MedLlama2+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_llama_32_3b"){
    fifth$method[iter]="LLama3.2_3B+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_phi4"){
    fifth$method[iter]="Phi-4+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_llama_33_70b"){
    fifth$method[iter]="Llama3.3_70B+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_MiniLM_L6_v2"){
    fifth$method[iter]="all-MiniLM-L6-v2+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_mpnet_base"){
    fifth$method[iter]="all-mpnet-base-v2+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_e5_large"){
    fifth$method[iter]="e5-large+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_gtr_t5_large"){
    fifth$method[iter]="gtr-t5-large+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_roberta"){
    fifth$method[iter]="all-roberta-large-v1+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_MiniLM_L12_v2"){
    fifth$method[iter]="all-MiniLM-L12-v2+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_labse"){
    fifth$method[iter]="LaBSE+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_sciBERT"){
    fifth$method[iter]="sciBERT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_sapBERT"){
    fifth$method[iter]="sapBERT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_cohere"){
    fifth$method[iter]="embed-english-v2.0+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_deepseek"){
    fifth$method[iter]="DeepSeek_8B+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_BioGPT"){
    fifth$method[iter]="BioGPT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_clincalBERT"){
    fifth$method[iter]="ClinicalBERT+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_e5large_v2"){
    fifth$method[iter]="e5large_v2+Euclidean Dist"
    
  }
  else if(fifth$method[iter]=="valid_euclidean_dist_nomic"){
    fifth$method[iter]="nomic+Euclidean Dist"
    
  }
}
rownames(all)<- NULL
rownames(fifth)<-NULL
print(all)
print(fifth)
