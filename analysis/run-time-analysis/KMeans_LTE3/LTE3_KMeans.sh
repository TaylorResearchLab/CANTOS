#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

Rscript 04B_5thed.R
Rscript 04B_all.R

Rscript 05B_5thed.R
Rscript 05B_all.R