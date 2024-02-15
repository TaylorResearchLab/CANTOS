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
