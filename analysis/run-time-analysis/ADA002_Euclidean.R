method_name <- "ADA002+Euclidean Distance"  
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


CT_embedding_df <- read.csv(paste(data_dir,"/CT_Embeddings_ADA2.csv",sep=""))
CT_embedding_df<-CT_embedding_df[,c(-1)]
# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
#WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <- read.csv(paste(data_dir,"/WHO_Aggregate_ADA2.csv",sep="")) #
WHO_embedding_df<-WHO_embedding_df[,c(-1)]

NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
#WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding


colnames(NCIT_embedding_df)[1]<-"Tumor_Names"
# Combined Embedding File for PC and Cluster Analysis
colnames(CT_embedding_df)<-colnames(NCIT_embedding_df) # CT embedding columns need to be fixed
colnames(WHO_embedding_df)<-colnames(NCIT_embedding_df)

combined_embedding_df <- rbind(CT_embedding_df,NCIT_embedding_df,WHO_embedding_df)
combined_embedding_df$Tumor_Names<- tolower(combined_embedding_df$Tumor_Names)
combined_embedding_df <- as.data.frame(combined_embedding_df %>% group_by(Tumor_Names) %>% summarise(across(everything(), list(mean))))
colnames(combined_embedding_df)<-colnames(NCIT_embedding_df)
rownames(combined_embedding_df)<-combined_embedding_df$Tumor_Names



distance_matrix_ADA002<- compute_embedding_distance(combined_embedding_df[1:nrow(combined_embedding_df),2:ncol(combined_embedding_df)],"euclidean")
rownames(distance_matrix_ADA002)<-combined_embedding_df$Tumor_Names
colnames(distance_matrix_ADA002)<-combined_embedding_df$Tumor_Names
end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time_2 <- Sys.time()
distance_matrix_ADA002<-as.data.frame(distance_matrix_ADA002) 
distance_matrix_ADA002_all <- distance_matrix_ADA002 %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_ADA002_5th <- distance_matrix_ADA002 %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_ADA002_ncit <- distance_matrix_ADA002 %>% dplyr::select(all_of(tolower(NCIT_embedding_df$Tumor_Names)))

ADA002_match_all<- nearest_match_embeddings(distance_matrix_ADA002_all,"ADA002") 
ADA002_match_5th<- nearest_match_embeddings(distance_matrix_ADA002_5th,"ADA002") 
ADA002_match_ncit<- nearest_match_embeddings(distance_matrix_ADA002_ncit,"ADA002") 
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



save.image(paste(analysis_dir,"/run-time-analysis/ada2_euclidean.RData",sep=""))
