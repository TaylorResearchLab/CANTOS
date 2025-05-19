#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

#Rscript lv_match.R
#Rscript cosine_match.R
#Rscript jw_match.R

#Rscript Medllama2_Euclidean.R
#Rscript llama32_Euclidean.R
#Rscript Phi4_Euclidean.R
#Rscript llama33_70b_Euclidean.R

#Rscript e5_large_Euclidean.R

#Rscript e5_large_v2_Euclidean.R
#Rscript cohere_euclidean.R

#Rscript nomic_Euclidean.R
#Rscript roberta_Euclidean.R
#Rscript modernBERT_Euclidean.R

#Rscript deepseek_8b_Euclidean.R
Rscript BioGPT_Euclidean.R
Rscript all-mpnet-base-v2_Euclidean.R
Rscript clinicalBERT_Euclidean.R
