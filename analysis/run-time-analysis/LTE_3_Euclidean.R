method_name <- "LTE-3+Euclidean Distance"  
start_time <- Sys.time()
# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(qdapRegex)
  library(jsonlite)
  library(httr)
  library(biomaRt)
  library(ghql)
  library(readxl)
  library(rollama)
  library(doParallel)
  library(foreach)
  library(stringdist)
  library(rdist)
  library(readr)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
source(paste(util_dir,"/compute_embedding_distance.R",sep=""))
source(paste(util_dir,"/nearest_match_embeddings.R",sep=""))


#
WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")


# Read CT embedding file 
ct_tumor_embeddings_v3_df<-read.csv(paste(data_dir,"/CT_Embeddings_V3.csv",sep=""))
ct_tumor_embeddings_v3_df<-ct_tumor_embeddings_v3_df[,c(-1)]

WHO_Terms_V3 <- read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep="")) #
WHO_Terms_V3<-WHO_Terms_V3[,c(-1)]

NCIT_embedding_v3<-read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep="")) 
NCIT_embedding_v3<-NCIT_embedding_v3[,c(-1)]

combined_embeddings_df<-rbind(ct_tumor_embeddings_v3_df,WHO_Terms_V3,NCIT_embedding_v3)#17660
combined_embeddings_df <- combined_embeddings_df%>%group_by(Tumor_Names) %>% summarise(across(everything(), list(mean)))#16721
combined_embeddings_df<-as.data.frame(combined_embeddings_df)
rownames(combined_embeddings_df)<-combined_embeddings_df$Tumor_Names
colnames(combined_embeddings_df)<-colnames(NCIT_embedding_v3)




distance_matrix_LTE_3<- compute_embedding_distance(combined_embeddings_df[1:nrow(combined_embeddings_df),2:ncol(combined_embeddings_df)],"euclidean")
rownames(distance_matrix_LTE_3)<-combined_embeddings_df$Tumor_Names
colnames(distance_matrix_LTE_3)<-combined_embeddings_df$Tumor_Names
end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time_2 <- Sys.time()
distance_matrix_LTE_3<-as.data.frame(distance_matrix_LTE_3) 
distance_matrix_LTE_3_all <- distance_matrix_LTE_3 %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_LTE_3_5th <- distance_matrix_LTE_3 %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_LTE_3_ncit <- distance_matrix_LTE_3 %>% dplyr::select(all_of(NCIT_embedding_v3$Tumor_Names))

LTE_3_match_all<- nearest_match_embeddings(distance_matrix_LTE_3_all,"lte_3") 
LTE_3_match_5th<- nearest_match_embeddings(distance_matrix_LTE_3_5th,"lte_3") 
LTE_3_match_ncit<- nearest_match_embeddings(distance_matrix_LTE_3_ncit,"lte_3") 
end_time_2 <- Sys.time()

elapsed_time_2<- as.numeric(difftime(end_time_2, start_time_2, units = "secs"))

time_elapsed_total<- elapsed_time_2+elapsed_time

benchmark <- data.frame(
  Method = method_name,
  Runtime_sec = round(time_elapsed_total, 2),
  stringsAsFactors = FALSE
)

output_file <- paste(analysis_dir,"/run-time-analysis/runtime_benchmarks_all_methods.csv",sep="")
if (!file.exists(output_file)) {
  write.csv(benchmark, output_file, row.names = FALSE)
} else {
  write.table(benchmark, output_file, append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}



save.image(paste(analysis_dir,"/run-time-analysis/lte-3_euclidean.RData",sep=""))
