# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
})


# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")


source(paste(util_dir,"/split_drugs.R",sep = ""))


# Load CT data
load(paste(input_dir,"/clinical_data.RData",sep=""))
load(paste(input_dir,"/eligibility_data.RData",sep=""))
load(paste(input_dir,"/conditions_data.RData",sep=""))
load(paste(input_dir,"/browse_conditions_data.RData",sep=""))

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


# CT-Disease-File
ct_disease_df <- as.data.frame(clinical_data$diseases)
colnames(ct_disease_df)<- "diseases"
ct_disease_df<-distinct(ct_disease_df)

# Write disease file 
write.csv(ct_disease_df,paste(intermediate_dir,"/ct_disease_df.csv",sep=""))