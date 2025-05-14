#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

#Rscript lv_match.R
#Rscript cosine_match.R
#Rscript jw_match.R

Rscript Levenshtein_AP.R