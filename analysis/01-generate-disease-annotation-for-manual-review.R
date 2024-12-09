# Annotate the 50K diseases automatically as cancer or not. 

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(tidyverse)
  
})


# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")


source(paste(util_dir,"/is_cancer_who.R",sep = ""))
source(paste(util_dir,"/create_disease_who_map.R",sep = ""))
source(paste(util_dir,"/cancer_fuzzy_match.R",sep = ""))
source(paste(util_dir,"/is_pedcan_who.R",sep = ""))
source(paste(util_dir,"/create_cancer_who_map.R",sep = ""))
source(paste(util_dir,"/pedcan_fuzzy_match.R",sep = ""))


# Read CT Diseases
ct_disease_df <- read.csv(paste(intermediate_dir,"/ct_disease_df.csv",sep=""))
ct_disease_df<-ct_disease_df %>% dplyr::select(nct_id,diseases)
  
who_general_words <- read_excel(paste(data_dir,"/who_cancer_key_words_general_5th_edition.xlsx",sep="")) # WHO Tumour List general (ADULT + PED)
who_ped_words<-  read_excel(paste(data_dir,"/who_cancer_key_words_paediatric_5th_edition.xlsx",sep="")) # WHO Tumour List PED

who_general_words$who_cancers<-tolower(who_general_words$who_cancers)
who_ped_words$who_paediatric_cancers<-tolower(who_ped_words$who_paediatric_cancers)


cancer_search_terms <- "cancer|carcinoma|adenocarcinoma|tumor|lymphoma|blast|myeloma|melanoma|leukemia|astrocytoma|malignant|neoplasm|neoplasia|
mesothelioma|ependymoma|glioma|thymoma|waldenstrom macroglobulinemia|myelodysplastic syndrome|polycythemia vera|myelofibrosis
|myeloproliferative|sarcoma|gist-plus syndrome|macroglobulinemia|mycosis fungoides|sezary's disease|plasmacytoma"



ind_cancer_search_term <- grepl(cancer_search_terms,ct_disease_df$diseases)

ct_disease_df$cancer_search_term<-NA

ct_disease_df$cancer_search_term[ind_cancer_search_term]<-"Yes"
ct_disease_df$cancer_search_term[!ind_cancer_search_term]<-"No"


#ct_disease_df<-is_cancer_who(ct_disease_df,who_general_words)
ct_disease_df<-is_cancer_who(ct_disease_df,who_general_words)
ct_disease_df<-ct_disease_df %>% dplyr::select(nct_id,diseases,cancer_search_term,is_cancer_who)

# Write this file and manually annotate with a column for validated_cancer_tumor ,PedCanTumor,Ped_Evidence 
write.csv(ct_disease_df,paste(intermediate_dir,"/ct_disease_df_manually_annotate.csv",sep="")) #16116 with either cancer_search_term =="Yes" or is_cancer_who=="Yes"

# load the  the manually annotated file
ct_disease_annot_adult_ped_df<-read.csv(paste(input_dir,"/tumor_annotated_no_ped.csv",sep=""))
ct_disease_annot_adult_ped_df<-ct_disease_annot_adult_ped_df[,c(-1)]

who_ped_words<-read_xlsx(paste(data_dir,"/WHO_Tumors/Paediatric_Classification_List.xlsx",sep=""))
colnames(who_ped_words)<-"who_paediatric_cancers"
who_ped_words$who_paediatric_cancers<-tolower(who_ped_words$who_paediatric_cancers)
ct_disease_annot_adult_ped_df<-is_pedcan_who(ct_disease_annot_adult_ped_df,who_ped_words)

write.csv(ct_disease_annot_adult_ped_df,paste(intermediate_dir,"/ct_disease_df_manually_annotate_pediatric.csv",sep="")) #16116 with either cancer_search_term =="Yes" or is_cancer_who=="Yes"
