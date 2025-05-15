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


tumor_all_gt<-read.csv(paste(result_dir,"/tumor_manually_validated_all.csv",sep = ""))
tumor_5thed_gt<-read.csv(paste(result_dir_5th,"/tumor_manually_validated_5th.csv",sep = ""))

tumor_all_gt<-tumor_all_gt[,c(-1)]
tumor_5thed_gt<-tumor_5thed_gt[,c(-1)]

tumor_all_gt<-tumor_all_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))
tumor_5thed_gt<-tumor_5thed_gt%>%filter(ground_truth %in% c("G","MG","G-Manual"))

correct_5thed_matrix<-tumor_5thed_gt[,seq(6,76,2)]
correct_all_matrix<-tumor_all_gt[,seq(6,76,2)]

colnames(correct_5thed_matrix)<-column_rename(colnames(correct_5thed_matrix))
colnames(correct_all_matrix)<-column_rename(colnames(correct_all_matrix))

result_5thed <- run_mcnemar_matrix(correct_5thed_matrix,adjust_method = "BH")
result_all <- run_mcnemar_matrix(correct_all_matrix,adjust_method = "BH")

adj_pvals_5thed<-result_5thed$adjusted_p_values
adj_pvals_all<- result_all$adjusted_p_values




signif_matrix_5thed <- ifelse(adj_pvals_5thed < 0.05, "*", "")
signif_matrix_all <- ifelse(adj_pvals_all < 0.05, "*", "")


#Heatmap for WHO-5th
pheatmap(adj_pvals_5thed,
         display_numbers = signif_matrix_5thed,
         main = "Pairwise McNemar Test: FDR-adjusted p-values , WHO-5th ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 10,
         legend = TRUE,
         na_col = "white") 

#Heatmap for WHO-all
pheatmap(adj_pvals_all,
         display_numbers = signif_matrix_all,
         main = "Pairwise McNemar Test: FDR-adjusted p-values , WHO-all ",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 10,
         legend = TRUE,
         na_col = "white") 


ci_results_5thed <- bootstrap_accuracy_ci_all_models(correct_5thed_matrix)
ci_results_all <- bootstrap_accuracy_ci_all_models(correct_all_matrix)

