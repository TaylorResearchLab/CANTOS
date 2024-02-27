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
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")


# Read CT embedding file 
embedding_df <- read.csv(paste(data_dir,"/disease_embeddings.csv",sep=""))
embedding_df<-embedding_df[order(embedding_df$Disease),]


# Read the annotated file
ct_disease_df <- read.csv(paste(input_dir,"/CT-Aug22-2023-Disease-File - clinical_trial_disease_aug_22_2023.csv",sep=""))
ct_tumor_df<- ct_disease_df %>% filter(validated_cancer_tumor=="Yes") %>%dplyr::select(diseases)



# Read missing CT embeddings 
CT_missing_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/missed_CT_result_embeddings.csv",sep=""))

# Combine all embeddings 

CT_embedding_df <- rbind(embedding_df,CT_missing_embedding_df)
CT_embedding_df <- unique(CT_embedding_df) # STILL CONTAINS DISEASES WITH EXACT STRING NAME BUT SLIGHTLY VARYING EMBEDDINGS
CT_embedding_df <- CT_embedding_df %>% group_by(Disease) %>% summarise(across(everything(), list(mean)))

# Match CT embeddings to CT files, some disease names in CT_embedding files are extra
#CT_only_df <- read.csv(paste(data_dir,"/CT_Only.csv",sep=""))
CT_embedding_agg_df <- ct_tumor_df %>% left_join(CT_embedding_df, by=c("diseases"="Disease"))
CT_embedding_agg_df<-CT_embedding_agg_df[order(CT_embedding_agg_df$diseases),] # Order by names so the rows with missing values are next to their matched rows 

# Some embeddings are not added due to name discrepecy,  commas were removed added when file was parsed in Open AI

ind_missing_commas <- which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)
ind_missing_commas <- unique(ind_missing_commas[,c(-2)])
ind_missing_commas_iteration_2 <-list()

for (iter in ind_missing_commas){
  print(iter)
  disease_name <- CT_embedding_agg_df$diseases[iter]
  disease_name_comma_removed <- gsub(",","",disease_name)
  ind_replacement <- which(CT_missing_embedding_df$Disease==disease_name_comma_removed)
  if(length(ind_replacement)>0){
    CT_embedding_agg_df[iter,2:1537] <- CT_missing_embedding_df[ind_replacement,2:1537]
  }
  else{
    disease_name_comma_removed <- gsub(","," ",disease_name)
    ind_replacement <- which(CT_missing_embedding_df$Disease==disease_name_comma_removed)
    if(length(ind_replacement)>0){
      CT_embedding_agg_df[iter,2:1537] <- CT_missing_embedding_df[ind_replacement,2:1537]
    }else{
      ind_missing_commas_iteration_2 <- c(ind_missing_commas_iteration_2,iter)
    }
    
  }
}


# Disease names with multiple commas need to be given embeddings


ind_missing_commas_iteration_2 <- unlist(ind_missing_commas_iteration_2)
tumor_list <- CT_missing_embedding_df$Disease
tumor_list <- gsub(",","",tumor_list)
tumor_list <- gsub(" ","",tumor_list)

for(iter in ind_missing_commas_iteration_2){
  disease_name <- CT_embedding_agg_df$diseases[iter]
  disease_name_comma_removed <- gsub(",","",disease_name)
  disease_name_comma_removed<- gsub(" ","",disease_name_comma_removed)
  ind_replacement <- which(tumor_list==disease_name_comma_removed)
  if(length(ind_replacement)>0){
    CT_embedding_agg_df[iter,2:1537]  <-CT_missing_embedding_df[ind_replacement,2:1537]
  }
}

# add embedding for 2 tumors with special characters in their names 
ind_missing_commas_iteration_3 <-  which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)
#ind_missing_commas_iteration_3 <- unique(ind_missing_commas_iteration_3[,1])
#CT_embedding_agg_df[ind_missing_commas_iteration_3[1],2:1537]<-CT_missing_embedding_df[1085,2:1537]
#CT_embedding_agg_df[ind_missing_commas_iteration_3[2],2:1537]<-CT_missing_embedding_df[1228,2:1537]


#Should print 0, then we are good to go
print(length(which(is.na(CT_embedding_agg_df[,2:1537]), arr.ind=TRUE)))

CT_embedding_agg_df<-unique(CT_embedding_agg_df) 

# Read NCIT Terms and WHO Terms with embedding and join them to the rest of the embedding list 
NCIT_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/NCIT_Neoplasm_Core_terms_text-embedding-ada-002_embeddings.csv",sep=""))
WHO_embedding_df <-read.csv(paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms_text-embedding-ada-002_embeddings.csv",sep=""))

NCIT_embedding_df<-NCIT_embedding_df[c(-1),] # Remove the header (column name) embedding
WHO_embedding_df<-WHO_embedding_df[c(-1),] # Remove the header (column name) embedding

rownames(NCIT_embedding_df)<-NULL
rownames(WHO_embedding_df)<-NULL

# Combined Embedding File for PC and Cluster Analysis
colnames(CT_embedding_agg_df)<-colnames(NCIT_embedding_df) # CT embedding columns need to be fixed
combined_embedding_df <- rbind(CT_embedding_agg_df,NCIT_embedding_df,WHO_embedding_df)
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

print(sum(ff[2,1:135]))

# Select the top 135 PCs
disease_transform = as.data.frame(-results$x[,1:135])

# Save this file "
write.csv(disease_transform,file =paste(intermediate_dir,"/disease_transform_pca.csv",sep="") )

