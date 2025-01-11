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
setwd(getwd())
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
util_dir <- file.path(root_dir, "util")
data_dir <- file.path(root_dir,"data")
input_dir <- file.path(root_dir,"data","input")
analysis_dir <- file.path(root_dir,"analysis")


# Load results to draw 1600 samples:
WHO_Results_all<- read.csv(paste(analysis_dir,"/results/WHO_Results_all_ose.csv",sep=""))
WHO_Results_5th<- read.csv(paste(analysis_dir,"/results_5th/WHO_Results_5thed_ose.csv",sep=""))

WHO_Results_all<-WHO_Results_all[,c(-1)]
WHO_Results_5th<-WHO_Results_5th[,c(-1)]

# Tumor sample
tumor_WHO_ALL <- sample_n(WHO_Results_all, 1600)
tumor_WHO_5th <- sample_n(WHO_Results_5th, 1600)

write.csv(tumor_WHO_ALL,paste(analysis_dir,"/results/tumor_all_validate_this_manually.csv",sep = ""))
write.csv(tumor_WHO_5th,paste(analysis_dir,"/results_5th/tumor_5th_validate_this_manually.csv",sep = ""))


