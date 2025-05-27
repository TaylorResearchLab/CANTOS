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
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

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




signif_matrix_5thed <- ifelse(adj_pvals_5thed < 0.05, "*", "")
signif_matrix_all <- ifelse(adj_pvals_all < 0.05, "*", "")


#Heatmap for WHO-5th
pheatmap(adj_pvals_5thed,
         display_numbers = signif_matrix_5thed,
         main = "FDR-Adjusted Pairwise McNemar’s p-values, WHO-5th (n=1044) ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 18,
         legend = TRUE,
         na_col = "white") 

#Heatmap for WHO-all
pheatmap(adj_pvals_all,
         display_numbers = signif_matrix_all,
         main = "FDR-Adjusted Pairwise McNemar’s p-values, WHO-all (n=1118) ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 18,
         legend = TRUE,
         na_col = "white") 


ci_results_5thed <- bootstrap_accuracy_ci_all_models(correct_5thed_matrix)
ci_results_all <- bootstrap_accuracy_ci_all_models(correct_all_matrix)



high_accuracy_methods<-c("LTE-3+AP",
                         "ADA-002+AP",
                         "LTE-3+KMeans",
                         "ADA-002+KMeans",
                         "LTE-3+Euclidean Distance",
                         "ADA-002+Euclidean Distance",
                         "all-MiniLM-L6-v2+Euclidean Distance",
                         "all-mpnet-base-v2+Euclidean Distance",
                         "all-MiniLM-L12-v2+Euclidean Distance",
                         "e5-large+Euclidean Distance",
                         "nomic+Euclidean Distance")

adj_all<-as.data.frame(adj_pvals_all)
adj_5th<-as.data.frame(adj_pvals_5thed)
adj_all_high<-adj_all%>%filter(rownames(adj_all)%in%high_accuracy_methods)
adj_all_high<-adj_all_high%>%select(any_of(high_accuracy_methods))

adj_5th_high<-adj_5th%>%filter(rownames(adj_5th)%in%high_accuracy_methods)
adj_5th_high<-adj_5th_high%>%select(any_of(high_accuracy_methods))


