#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 


Rscript 04A-ADA002_AP_5th.R
Rscript 04A-ADA002_all.R

Rscript 06A_ADA002_5th.R
Rscript 06A_ADA002_all.R

Rscript 07A_ADA002_5th.R
Rscript 07A_ADA002_all.R

Rscript 08_ADA002_5th.R
Rscript 08_ADA002_all.R

Rscript 09_ADA002_5th.R
Rscript 09_ADA002_all.R

