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
  library(pheatmap)
  library(grid)
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

source(paste(util_dir,"/run_mcnemar_matrix.R",sep=""))
source(paste(util_dir,"/column_rename.R",sep=""))
source(paste(util_dir,"/bootstrap_accuracy_ci_all_models.R",sep=""))



tumor_all_gt<-read_excel(paste(result_dir,"/tumor_manually_validated_all_corrected_May23.xlsx",sep = ""))
tumor_5thed_gt<-read_excel(paste(result_dir_5th,"/tumor_manually_validated_5th_corrected_May23.xlsx",sep = ""))

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

correct_5thed_matrix<-tumor_5thed_gt[,seq(6,76,2)]
correct_all_matrix<-tumor_all_gt[,seq(6,76,2)]

colnames(correct_5thed_matrix)<-column_rename(colnames(correct_5thed_matrix))
colnames(correct_all_matrix)<-column_rename(colnames(correct_all_matrix))

correct_5thed_matrix<-as.matrix(correct_5thed_matrix)
correct_all_matrix<-as.matrix(correct_all_matrix)

result_5thed <- run_mcnemar_matrix(correct_5thed_matrix,adjust_method = "BH")
result_all <- run_mcnemar_matrix(correct_all_matrix,adjust_method = "BH")

adj_pvals_5thed<-result_5thed$adjusted_p_values
adj_pvals_all<- result_all$adjusted_p_values

rownames(adj_pvals_5thed)[which(rownames(adj_pvals_5thed)=="LTE-3+KMeans")]<-"LTE-3+K-means"
colnames(adj_pvals_5thed)[which(colnames(adj_pvals_5thed)=="LTE-3+KMeans")]<-"LTE-3+K-means"

rownames(adj_pvals_5thed)[which(rownames(adj_pvals_5thed)=="ADA-002+KMeans")]<-"ADA-002+K-means"
colnames(adj_pvals_5thed)[which(colnames(adj_pvals_5thed)=="ADA-002+KMeans")]<-"ADA-002+K-means"


rownames(adj_pvals_all)[which(rownames(adj_pvals_all)=="LTE-3+KMeans")]<-"LTE-3+K-means"
colnames(adj_pvals_all)[which(colnames(adj_pvals_all)=="LTE-3+KMeans")]<-"LTE-3+K-means"

rownames(adj_pvals_all)[which(rownames(adj_pvals_all)=="ADA-002+KMeans")]<-"ADA-002+K-means"
colnames(adj_pvals_all)[which(colnames(adj_pvals_all)=="ADA-002+KMeans")]<-"ADA-002+K-means"



signif_matrix_5thed <- ifelse(adj_pvals_5thed < 0.05, "*", "")
signif_matrix_all <- ifelse(adj_pvals_all < 0.05, "*", "")


method_ranking_5th <- c("all-MiniLM-L12-v2+Euclidean distance", "e5-large+Euclidean distance", "all-mpnet-base-v2+Euclidean distance", "LTE-3+Euclidean distance", "all-MiniLM-L6-v2+Euclidean distance",
                        "LTE-3+AP", "LTE-3+K-means", "ADA-002+Euclidean distance", "ADA-002+AP", "ADA-002+K-means", 
                        "nomic+Euclidean distance", "all-roberta-large-v1+Euclidean distance", "e5large_v2+Euclidean distance",
                        "gtr-t5-large+Euclidean distance", "embed-english-v2.0+Euclidean distance", "sapBERT+Euclidean distance",
                        "BioGPT+Euclidean distance", "ClinicalBERT+Euclidean distance", "LaBSE+Euclidean distance",
                        "MedLlama_13B+Euclidean distance", "Levenshtein", "Levenshtein+AP", "Llama3.3_70B+Euclidean distance", 
                        "Jaro Winkler", "Jaro Winkler+AP", "Cosine", "Cosine+AP", "DeepSeek_8B+Euclidean distance", "PubMedBERT+Euclidean distance",
                        "LLama3.2_3B+Euclidean distance", "BioBERT+Euclidean distance", "LLama3.0+Euclidean distance",
                        "sciBERT+Euclidean distance", "ModernBERT+Euclidean distance", "Phi-4+Euclidean distance", 
                        "MedLlama2+Euclidean distance")





method_ranking_all<-c("LTE-3+Euclidean distance", "e5-large+Euclidean distance", "LTE-3+AP", "all-mpnet-base-v2+Euclidean distance", "ADA-002+Euclidean distance",
                      "all-MiniLM-L12-v2+Euclidean distance", "all-MiniLM-L6-v2+Euclidean distance", "ADA-002+AP", 
                      "LTE-3+K-means", "ADA-002+K-means", "nomic+Euclidean distance", "all-roberta-large-v1+Euclidean distance",
                      "gtr-t5-large+Euclidean distance", "e5large_v2+Euclidean distance", "embed-english-v2.0+Euclidean distance",
                      "sapBERT+Euclidean distance", "ClinicalBERT+Euclidean distance", "BioGPT+Euclidean distance", 
                      "LaBSE+Euclidean distance", "Levenshtein", "MedLlama_13B+Euclidean distance", "Levenshtein+AP",
                      "Llama3.3_70B+Euclidean distance", "Jaro Winkler", "Jaro Winkler+AP", "Cosine", 
                      "DeepSeek_8B+Euclidean distance", "Cosine+AP", "PubMedBERT+Euclidean distance", 
                      "sciBERT+Euclidean distance", "LLama3.2_3B+Euclidean distance", "BioBERT+Euclidean distance",
                      "LLama3.0+Euclidean distance", "ModernBERT+Euclidean distance", "Phi-4+Euclidean distance", 
                      "MedLlama2+Euclidean distance")


print(length(!which(!rownames(adj_pvals_5thed)%in%method_ranking_5th)))
print(length(!which(!rownames(adj_pvals_all)%in%method_ranking_all)))

ordered_adj_pval_5th <- adj_pvals_5thed[method_ranking_5th, method_ranking_5th]
ordered_signif_5th <- signif_matrix_5thed[method_ranking_5th, method_ranking_5th]

#Heatmap for WHO-5th
jpeg(paste(root_dir,"/Paper/MLWA/High_Res_Fig/figure8.jpg",sep=""), width = 40, height = 35, units = "cm", res = 300)
pheatmap(ordered_adj_pval_5th,
         display_numbers = ordered_signif_5th,
         main = "FDR-Adjusted Pairwise McNemar’s p-values, WHO-5th (n=1044) ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 16,            # overall font size
         fontsize_row = 16,        # row names (y-axis)
         fontsize_col = 16,        # column names (x-axis)
         fontsize_number = 16,     # text inside cells
         legend = TRUE,
         na_col = "white") 
dev.off()

ordered_adj_pval_all <- adj_pvals_all[method_ranking_all, method_ranking_all]
ordered_signif_all <- signif_matrix_all[method_ranking_all, method_ranking_all]



#Heatmap for WHO-all
jpeg(paste(root_dir,"/Paper/MLWA/High_Res_Fig/figure7.jpg",sep=""), width = 40, height = 35, units = "cm", res = 300)

pheatmap(ordered_adj_pval_all,
         display_numbers = ordered_signif_all,
         main = "FDR-Adjusted Pairwise McNemar’s p-values, WHO-all (n=1118) ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize = 16,            # overall font size
         fontsize_row = 16,        # row names (y-axis)
         fontsize_col = 16,        # column names (x-axis)
         fontsize_number = 16,     # text inside cells         legend = TRUE,
         na_col = "white") 


dev.off()



ci_results_5thed <- bootstrap_accuracy_ci_all_models(correct_5thed_matrix)
ci_results_all <- bootstrap_accuracy_ci_all_models(correct_all_matrix)



high_accuracy_methods<-c("LTE-3+AP",
                         "ADA-002+AP",
                         "LTE-3+K-means",
                         "ADA-002+K-means",
                         "LTE-3+Euclidean distance",
                         "ADA-002+Euclidean distance",
                         "all-MiniLM-L6-v2+Euclidean distance",
                         "all-mpnet-base-v2+Euclidean distance",
                         "all-MiniLM-L12-v2+Euclidean distance",
                         "e5-large+Euclidean distance",
                         "nomic+Euclidean distance")

adj_all<-as.data.frame(ordered_adj_pval_all)
adj_5th<-as.data.frame(ordered_adj_pval_5th)
adj_all_high<-adj_all%>%filter(rownames(adj_all)%in%high_accuracy_methods)
adj_all_high<-adj_all_high%>%select(any_of(high_accuracy_methods))

adj_5th_high<-adj_5th%>%filter(rownames(adj_5th)%in%high_accuracy_methods)
adj_5th_high<-adj_5th_high%>%select(any_of(high_accuracy_methods))

colnames(adj_all_high)[(which(adj_all_high[5,]>0.05))]


biomed_methods<-c("sapBERT+Euclidean distance", "ClinicalBERT+Euclidean distance", "BioGPT+Euclidean distance",
                  "MedLlama_13B+Euclidean distance", "PubMedBERT+Euclidean distance", "sciBERT+Euclidean distance", 
                  "BioBERT+Euclidean distance",  "MedLlama2+Euclidean distance")

adj_all_biomed<-adj_all%>%filter(rownames(adj_all)%in%c(biomed_methods,high_accuracy_methods))
adj_all_biomed<-adj_all_biomed%>%select(any_of(c(biomed_methods,high_accuracy_methods)))


adj_5th_biomed<-adj_5th%>%filter(rownames(adj_5th)%in%c(biomed_methods,high_accuracy_methods))
adj_5th_biomed<-adj_5th_biomed%>%select(any_of(c(biomed_methods,high_accuracy_methods)))

colnames(adj_5th_high[which(adj_5th_high[9,]>=0.05)])
colnames(adj_5th_high[which(adj_5th_high[9,]<0.05)])

write.csv(ordered_adj_pval_5th,paste(result_dir_5th,"/mcnemars_adj_pval_5th.csv",sep=""))
write.csv(ordered_adj_pval_all,paste(result_dir,"/mcnemars_adj_pval_all.csv",sep=""))
