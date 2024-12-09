suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(tidyverse)
  library(stringi)
  library(readxl)
  library(openxlsx)
})

# Set the directories
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"input")
analysis_dir <- file.path(root_dir,"analysis")
intermediate_dir <- file.path(analysis_dir,"intermediate")
result_dir <- file.path(analysis_dir,"results")

# Load WHO Data 3rd edition
excel_sheet_3rd_edition <- loadWorkbook((paste(data_dir,"/WHO_Tumors/3rd_edition/WHO_3rd_edition.xlsx",sep="")))
sheetNames_3rd_edition <- sheets(excel_sheet_3rd_edition)
df_3rd_edition<-data.frame()

for(i in 1:length(sheetNames_3rd_edition))
{
  df_3rd_edition<- rbind(readWorkbook(excel_sheet_3rd_edition,sheet = i,colNames = FALSE),df_3rd_edition)
}

# Load WHO Data 4th edition
excel_sheet_4th_edition <- loadWorkbook((paste(data_dir,"/WHO_Tumors/4th_edition/WHO_Tumor_4th_edition.xlsx",sep="")))
sheetNames_4th_edition <- sheets(excel_sheet_4th_edition)

df_4th_edition<-data.frame()
for(i in 1:length(sheetNames_4th_edition))
{
  df_4th_edition<- rbind(readWorkbook(excel_sheet_4th_edition,sheet = i,colNames = FALSE),df_4th_edition)
}

# Load WHO Data 5th edition
df_5th_edition <- read.csv((paste(data_dir,"/dt_input_file_6_dec/WHO_Only_terms.csv",sep="")))


#rename the columns
colnames(df_3rd_edition)<-"Tumor_Names"
colnames(df_4th_edition)<-"Tumor_Names"
colnames(df_5th_edition)<-"Tumor_Names"

#lower case , and trim white spaces to left
df_3rd_edition$Tumor_Names<- tolower(df_3rd_edition$Tumor_Names)
df_4th_edition$Tumor_Names<- tolower(df_4th_edition$Tumor_Names)
df_5th_edition$Tumor_Names<- tolower(df_5th_edition$Tumor_Names)

df_3rd_edition$Tumor_Names<- str_trim(df_3rd_edition$Tumor_Names)
df_4th_edition$Tumor_Names<- str_trim(df_4th_edition$Tumor_Names)
df_5th_edition$Tumor_Names<- str_trim(df_5th_edition$Tumor_Names)

# distinct
df_3rd_edition<- distinct(df_3rd_edition)
df_4th_edition<- distinct(df_4th_edition)
df_5th_edition<- distinct(df_5th_edition)


# and, /, "," in the text 
df_3rd_edition <- df_3rd_edition %>% mutate(is_AND= case_when(str_detect(Tumor_Names,"and")~"Yes",TRUE~"No"))
df_3rd_edition <- df_3rd_edition %>% mutate(is_slash= case_when(str_detect(Tumor_Names,"/")~"Yes",TRUE~"No"))
df_3rd_edition <- df_3rd_edition %>% mutate(is_comma= case_when(str_detect(Tumor_Names,",")~"Yes",TRUE~"No"))
df_3rd_edition<-df_3rd_edition %>% arrange(is_AND, is_slash,is_comma) 

# and, /, "," in the text 
df_4th_edition <- df_4th_edition %>% mutate(is_AND= case_when(str_detect(Tumor_Names,"and")~"Yes",TRUE~"No"))
df_4th_edition <- df_4th_edition %>% mutate(is_slash= case_when(str_detect(Tumor_Names,"/")~"Yes",TRUE~"No"))
df_4th_edition <- df_4th_edition %>% mutate(is_comma= case_when(str_detect(Tumor_Names,",")~"Yes",TRUE~"No"))
df_4th_edition<-df_4th_edition %>% arrange(is_AND, is_slash,is_comma) 

# and, /, "," in the text 
df_5th_edition <- df_5th_edition %>% mutate(is_AND= case_when(str_detect(Tumor_Names,"and")~"Yes",TRUE~"No"))
df_5th_edition <- df_5th_edition %>% mutate(is_slash= case_when(str_detect(Tumor_Names,"/")~"Yes",TRUE~"No"))
df_5th_edition <- df_5th_edition %>% mutate(is_comma= case_when(str_detect(Tumor_Names,",")~"Yes",TRUE~"No"))
df_5th_edition<-df_5th_edition %>% arrange(is_AND, is_slash,is_comma) 
# Write the csv
write.xlsx(df_5th_edition,paste(data_dir,"/WHO_Tumors/intermediate/df_5th_edition_manual_edit.xlsx",sep=""))
write.xlsx(df_4th_edition,paste(data_dir,"/WHO_Tumors/intermediate/df_4th_edition_manual_edit.xlsx",sep=""))
write.xlsx(df_3rd_edition,paste(data_dir,"/WHO_Tumors/intermediate/df_3rd_edition_manual_edit.xlsx",sep=""))



#### Read post edit files
WHO_df_3rd<- read.xlsx(paste(data_dir,"/WHO_Tumors/intermediate/df_3rd_edition_post_edit.xlsx",sep=""),colNames = FALSE)
WHO_df_4th<- read.xlsx(paste(data_dir,"/WHO_Tumors/intermediate/df_4th_edition_post_edit.xlsx",sep=""),colNames = FALSE)
WHO_df_5th<- read.xlsx(paste(data_dir,"/WHO_Tumors/intermediate/df_5th_edition_post_edit.xlsx",sep=""),colNames = FALSE)

#lower case , and trim white spaces to left
colnames(WHO_df_3rd)<- "Tumor_Names"
colnames(WHO_df_4th)<- "Tumor_Names"
colnames(WHO_df_5th)<- "Tumor_Names"

WHO_df_3rd$Tumor_Names<- tolower(WHO_df_3rd$Tumor_Names)
WHO_df_4th$Tumor_Names<- tolower(WHO_df_4th$Tumor_Names)
WHO_df_5th$Tumor_Names<- tolower(WHO_df_5th$Tumor_Names)

WHO_df_3rd$Tumor_Names<- str_trim(WHO_df_3rd$Tumor_Names)
WHO_df_4th$Tumor_Names<- str_trim(WHO_df_4th$Tumor_Names)
WHO_df_5th$Tumor_Names<- str_trim(WHO_df_5th$Tumor_Names)

# all who tumors combined

WHO_Tumors <- rbind(WHO_df_3rd,WHO_df_4th,WHO_df_5th)
WHO_Tumors<- distinct(WHO_Tumors)


# Add editions to who tumors
WHO_Tumors <- WHO_Tumors %>% mutate(edition_5th =case_when(Tumor_Names %in% WHO_df_5th$Tumor_Names~"Yes", TRUE ~ "No"))
WHO_Tumors <- WHO_Tumors %>% mutate(edition_4th =case_when(Tumor_Names %in% WHO_df_4th$Tumor_Names~"Yes", TRUE ~ "No"))
WHO_Tumors <- WHO_Tumors %>% mutate(edition_3rd =case_when(Tumor_Names %in% WHO_df_3rd$Tumor_Names~"Yes", TRUE ~ "No"))

# Write file for WHO Tumor
write.xlsx(WHO_Tumors,paste(data_dir,"/WHO_Tumors/result/WHO_Tumor_all_edition.xlsx",sep=""))

WHO_OLD<- WHO_Tumors[which(!(WHO_Tumors$Tumor_Names %in% WHO_embedding_df$Disease)),]
write.csv(WHO_OLD,paste(data_dir,"/WHO_Tumors/result/WHO_OLD.csv",sep=""))