method_name <- "pubmedbert+Euclidean Distance"
start_time <- Sys.time()

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(readxl)
  library(rdist)
  library(readr)
})
setwd(getwd())
root_dir<-"/home/lahiria/CANTOS-RUN-TIME"
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"data","input")
analysis_dir <- file.path(root_dir,"analysis")

source(paste(util_dir,"/compute_embedding_distance.R",sep=""))
source(paste(util_dir,"/nearest_match_embeddings.R",sep=""))



#
ct_terms <- read.csv(paste(data_dir,"/CT_Embeddings_V3.csv",sep = ""))
ncit_terms <- read.csv(paste(data_dir,"/NCIT_Embeddings_V3.csv",sep = ""))
who_terms <- read.csv(paste(data_dir,"/WHO_Terms_All_V3.csv",sep = ""))


ct_terms<- ct_terms %>% dplyr::select(Tumor_Names)
ncit_terms<- ncit_terms %>% dplyr::select(Tumor_Names)
who_terms<- who_terms %>% dplyr::select(Tumor_Names)

ct_terms2<-ct_terms %>% dplyr::pull(Tumor_Names)
who_terms2<-who_terms %>% dplyr::pull(Tumor_Names)
ncit_terms2<-ncit_terms %>% dplyr::pull(Tumor_Names)


tumor_terms<-c(ct_terms2,who_terms2,ncit_terms2)

print("pubmedbert Loading")
pubmedbert_embeddings <-read.csv(paste(data_dir,"/Embeddings/pubmedbert-base-embeddings.csv",sep="")) #768 dimension
end_time <- Sys.time()
elapsed_time<- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time<-Sys.time()
pubmedbert_embeddings<-pubmedbert_embeddings %>% group_by(Tumor_Names) %>% summarise_all("mean")
end_time <- Sys.time()
elapsed_time2<- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time<-Sys.time()

WHO_Terms_All <-readxl::read_xlsx(paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))
WHO_Terms_5th<-WHO_Terms_All%>%filter(edition_5th=="Yes")

print("pubmedbert Distance Matrix")
distance_matrix_pubmedbert<-compute_embedding_distance(pubmedbert_embeddings[1:nrow(pubmedbert_embeddings),2:ncol(pubmedbert_embeddings)],"euclidean")
rownames(distance_matrix_pubmedbert)<-pubmedbert_embeddings$Tumor_Names
colnames(distance_matrix_pubmedbert)<-pubmedbert_embeddings$Tumor_Names




end_time <- Sys.time()
elapsed_time3<- as.numeric(difftime(end_time, start_time, units = "secs"))

print("pubmedbert Distance Matrix Computed")

start_time<-Sys.time()
save.image("pubmedbert.RData")

distance_matrix_pubmedbert<-as.data.frame(distance_matrix_pubmedbert) # EXECUTED


distance_matrix_pubmedbert_all <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_All$Tumor_Names))
distance_matrix_pubmedbert_5th <- distance_matrix_pubmedbert %>% dplyr::select(all_of(WHO_Terms_5th$Tumor_Names))
distance_matrix_pubmedbert_ncit <- distance_matrix_pubmedbert %>% dplyr::select(all_of(ncit_terms$Tumor_Names))


pubmedbert_all<- nearest_match_embeddings(distance_matrix_pubmedbert_all,"pubmedbert") # Executed
pubmedbert_5th<- nearest_match_embeddings(distance_matrix_pubmedbert_5th,"pubmedbert") # Executed
pubmedbert_ncit<- nearest_match_embeddings(distance_matrix_pubmedbert_ncit,"pubmedbert") # Executed
end_time <- Sys.time()
elapsed_time4<- as.numeric(difftime(end_time, start_time, units = "secs"))

print("pubmedbert Computed")



objs <- ls(envir = .GlobalEnv)

# Get the size of each object
sizes <- sapply(objs, function(x) object.size(get(x, envir = .GlobalEnv)))

# Convert sizes to readable format and summarize
sizes_df <- data.frame(
  Object = objs,
  Size_MB = round(sizes / (1024^2), 3)
)

# Sort by size descending
sizes_df <- sizes_df[order(-sizes_df$Size_MB), ]

# Print the total memory used
cat("Total memory used in Global Environment:", round(sum(sizes) / (1024^2), 3), "MB\n")


time_elapsed_total<- elapsed_time+elapsed_time2+elapsed_time3+elapsed_time4

benchmark <- data.frame(
  Method = method_name,
  Runtime_sec = round(time_elapsed_total, 2),
  Memory_gb=round(sum(sizes) / (1024^2), 3)/1000,
  stringsAsFactors = FALSE
)

print(benchmark)

runtime_file <- read.csv(paste(analysis_dir,"/run-time-analysis/runtime_benchmarks_all_methods.csv",sep=""))
runtime_file<-runtime_file[,-1]
runtime_file <- rbind(runtime_file, benchmark)
write.csv(runtime_file, paste(analysis_dir,"/run-time-analysis/runtime_benchmarks_all_methods.csv",sep=""))

