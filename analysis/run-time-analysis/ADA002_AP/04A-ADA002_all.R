# Loads embeddings ADA2 for CT, WHO, NCIT  Tumors and then performs PCA. 
start_time_4a_all<-Sys.time()
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  #  library(qdapRegex)
  library(jsonlite)
  library(httr)
  # library(biomaRt)
  # library(ghql)
  library(readxl)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
run_time_analysis<-file.path(analysis_dir,"run-time-analysis")
intermediate_dir <- file.path(run_time_analysis,"intermediate")

# Read CT embedding file 
CT_embedding_df <- read.csv(paste(data_dir,"/CT_Embeddings_ADA2.csv",sep=""))
CT_embedding_df<-CT_embedding_df[,c(-1)]
# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
#WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <- read.csv(paste(data_dir,"/WHO_Aggregate_ADA2.csv",sep="")) #
WHO_embedding_df<-WHO_embedding_df[,c(-1)]

NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
#WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

# Combined Embedding File for PC and Cluster Analysis
colnames(CT_embedding_df)<-colnames(NCIT_embedding_df) # CT embedding columns need to be fixed
colnames(WHO_embedding_df)<-colnames(NCIT_embedding_df)

combined_embedding_df <- rbind(CT_embedding_df,NCIT_embedding_df,WHO_embedding_df)
combined_embedding_df$Disease<- tolower(combined_embedding_df$Disease)
combined_embedding_df <- as.data.frame(combined_embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean))))
colnames(combined_embedding_df)<-colnames(NCIT_embedding_df)
rownames(combined_embedding_df)<-combined_embedding_df$Disease
combined_embedding_df<-combined_embedding_df[,c(-1)]

# Now perform PCA 
results <- prcomp(combined_embedding_df, scale = TRUE)
eigs <- results$sdev^2
ff<-rbind(SD = sqrt(eigs),Proportion = eigs/sum(eigs), Cumulative = cumsum(eigs)/sum(eigs))

# Check the number of components needed to capture 80% variance at least

print(sum(ff[2,1:141]))
# Select the top 141 PCs
disease_transform = as.data.frame(-results$x[,1:141])

end_time_4a_all<-Sys.time()
elapsed_time_4a_all<- as.numeric(difftime(end_time_4a_all, start_time_4a_all, units = "secs"))

# Save this file "
write.csv(disease_transform,file =paste(intermediate_dir,"/disease_transform_pca_ada2.csv",sep="") )
save(combined_embedding_df,file = paste(intermediate_dir,"/combined_embedding_ada2_df.RData",sep=""))
#save.image(file = "script4a.RData")

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

Runtime_4A_all = round(elapsed_time_4a_all, 2)
Memory_gb_4a_all=round(sum(sizes) / (1024^2), 3)/1000

save.image(paste(analysis_dir,"/run-time-analysis/ADA002_AP/4A-ADA002_all.RData",sep=""))

