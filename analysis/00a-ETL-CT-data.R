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
  # For embedding
  library(word2vec)
  library(keras)
  library(reticulate)
  library(purrr)
  library(stringdist)
})


# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")

# Load functions
source(paste(util_dir,"/split_drugs.R",sep = ""))
source(paste(util_dir,"/create_disease_who_map.R",sep = ""))
source(paste(util_dir,"/clustered_tumors.R",sep = ""))
source(paste(util_dir,"/contains_target_word.R",sep = ""))
source(paste(util_dir,"/string_dissimilarity.R",sep = ""))

# Load CT data
load(paste(input_dir,"/clinical_data.RData",sep=""))
load(paste(input_dir,"/eligibility_data.RData",sep=""))
load(paste(input_dir,"/conditions_data.RData",sep=""))
load(paste(input_dir,"/browse_conditions_data.RData",sep=""))

df_tumor<- read.csv(paste(input_dir,"/ct_master_tumor_ped_annotated.csv",sep=""))


# Relationship between Drug and Disease in CT and NCT IDS
clinical_data <- clinical_data %>% dplyr::left_join(eligibility_data,by="nct_id")
clinical_data<-clinical_data %>% dplyr::select(nct_id,intervention_type,name)
clinical_data<-clinical_data %>% dplyr::filter(intervention_type %in% c("Drug","Biological","Combination Product","Genetic" ))
colnames(clinical_data)[3]<-"Intervention_CT"

clinical_data<-split_drugs(clinical_data)
clinical_data<-distinct(clinical_data)
clinical_data <- clinical_data %>% dplyr::select(nct_id,Intervention_CT)

conditions_data <- conditions_data %>% dplyr::select(nct_id,downcase_name)
browse_conditions_data<-browse_conditions_data %>% dplyr::select(nct_id,downcase_mesh_term)
clinical_data <- clinical_data %>% dplyr::left_join(conditions_data,by="nct_id") # contains CT drugs and disease info with NCT_IDs
colnames(clinical_data)[3]<-"diseases"

df_tumor_drug<- dplyr::left_join(df_tumor,clinical_data, by="diseases")
df_tumor_drug <- df_tumor_drug %>% dplyr::filter(validated_cancer_tumor!="No")
df_tumor_drug <- df_tumor_drug %>% dplyr::filter(disease_name_inconclusive!="Yes")

